from pathlib import Path
from collections import defaultdict, Counter
import csv
import re

CDS_BASE = Path("data/cds")
WEIGHT_DIR = Path("work/codon/species_weights")

THRESHOLD_FILE = Path(
    "work/codon/frozen_thresholds/species_low_adaptation_thresholds.tsv"
)

FAMILY_TABLE = Path(
    "work/topology/merged/topology_qualified_families.tsv"
)

PROTEIN_FASTA = Path(
    "work/topology/bacillales_ge8of10_proteins.faa"
)

OUTDIR = Path("work/codon/segments")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_SEGMENTS = OUTDIR / "low_adaptation_segments.tsv"
OUT_PROTEINS = OUTDIR / "protein_low_adaptation_summary.tsv"
OUT_SPECIES = OUTDIR / "species_low_adaptation_summary.tsv"

MIN_RUN = 3


def read_fasta(path):
    records = {}
    header = None
    seq = []

    with path.open() as f:
        for line in f:
            line = line.rstrip()

            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq)

                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())

        if header is not None:
            records[header] = "".join(seq)

    return records


# --------------------------------------------------
# 1. Frozen thresholds
# --------------------------------------------------

thresholds = {}

with THRESHOLD_FILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        thresholds[row["accession"]] = float(row["threshold"])

print("Frozen species thresholds:", len(thresholds))


# --------------------------------------------------
# 2. Qualified orthogroups
# --------------------------------------------------

qualified = set()

with FAMILY_TABLE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        if row["topology_qualified"] == "1":
            qualified.add(row["orthogroup"])

print("Qualified families:", len(qualified))


# --------------------------------------------------
# 3. Target proteins
# --------------------------------------------------

targets = {}

protein_records = read_fasta(PROTEIN_FASTA)

for header, protein_seq in protein_records.items():

    seq_id = header.split()[0]
    parts = seq_id.split("|", 2)

    if len(parts) != 3:
        continue

    og, accession, protein_id = parts

    if og not in qualified:
        continue

    targets[seq_id] = {
        "orthogroup": og,
        "accession": accession,
        "protein_id": protein_id,
        "protein_length": len(protein_seq),
    }

print("Target proteins:", len(targets))


# --------------------------------------------------
# 4. Load weights
# --------------------------------------------------

weights_by_species = {}

for accession in thresholds:

    wf = WEIGHT_DIR / f"{accession}_codon_weights.tsv"

    weights = {}

    with wf.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            weights[row["codon"]] = float(
                row["relative_weight"]
            )

    weights_by_species[accession] = weights


# --------------------------------------------------
# 5. Load CDS indexes
# --------------------------------------------------

cds_by_species = {}

for accession in thresholds:

    files = list(
        (CDS_BASE / accession).rglob(
            "cds_from_genomic.fna"
        )
    )

    if len(files) != 1:
        raise RuntimeError(
            f"{accession}: expected exactly 1 CDS file, "
            f"found {len(files)}"
        )

    protein_to_cds = {}

    for header, seq in read_fasta(files[0]).items():

        m = re.search(
            r"\[protein_id=([^\]]+)\]",
            header
        )

        if m:
            protein_to_cds[m.group(1)] = seq.upper()

    cds_by_species[accession] = protein_to_cds


# --------------------------------------------------
# 6. Call segments
# --------------------------------------------------

segment_rows = []
protein_rows = []

species_stats = defaultdict(lambda: {
    "proteins": 0,
    "proteins_with_segment": 0,
    "eligible_positions": 0,
    "low_positions": 0,
    "segments": 0,
})

for seq_id, rec in targets.items():

    og = rec["orthogroup"]
    accession = rec["accession"]
    protein_id = rec["protein_id"]
    protein_length = rec["protein_length"]

    cds = cds_by_species[accession].get(protein_id)

    if cds is None:
        raise RuntimeError(
            f"CDS missing: {seq_id}"
        )

    codons = [
        cds[i:i+3]
        for i in range(0, len(cds), 3)
    ]

    # Remove terminal stop from positional analysis,
    # but retain original 1-based protein/codon coordinates.
    if codons and codons[-1] in {"TAA", "TAG", "TGA"}:
        analysis_last = len(codons) - 1
    else:
        analysis_last = len(codons)

    threshold = thresholds[accession]
    weights = weights_by_species[accession]

    # Each item:
    # (original 1-based codon position, codon, weight, is_low)
    positions = []

    # Start from index 1 in Python = biological codon position 2.
    # Thus annotated initiation codon is omitted.
    for i in range(1, analysis_last):

        codon = codons[i]
        biological_pos = i + 1

        if codon not in weights:
            raise RuntimeError(
                f"{seq_id}: codon {codon} has no weight"
            )

        w = weights[codon]
        is_low = w <= threshold

        positions.append(
            (biological_pos, codon, w, is_low)
        )

    low_count = sum(x[3] for x in positions)

    # Detect consecutive low-adaptation runs
    runs = []
    current = []

    for item in positions:

        pos, codon, w, is_low = item

        if is_low:
            if current:
                prev_pos = current[-1][0]

                # Require true consecutive codon positions
                if pos != prev_pos + 1:
                    if len(current) >= MIN_RUN:
                        runs.append(current)
                    current = []

            current.append(item)

        else:
            if len(current) >= MIN_RUN:
                runs.append(current)

            current = []

    if len(current) >= MIN_RUN:
        runs.append(current)

    for segment_index, run in enumerate(runs, start=1):

        start = run[0][0]
        end = run[-1][0]
        center = (start + end) / 2

        segment_rows.append({
            "orthogroup": og,
            "accession": accession,
            "protein_id": protein_id,
            "sequence_id": seq_id,
            "segment_index": segment_index,
            "segment_start": start,
            "segment_end": end,
            "segment_center": center,
            "segment_length_codons": len(run),
            "codons": ";".join(x[1] for x in run),
            "weights": ";".join(
                f"{x[2]:.8f}" for x in run
            ),
            "species_threshold": threshold,
        })

    protein_rows.append({
        "orthogroup": og,
        "accession": accession,
        "protein_id": protein_id,
        "sequence_id": seq_id,
        "protein_length": protein_length,
        "eligible_positions": len(positions),
        "low_adaptation_positions": low_count,
        "low_adaptation_fraction": (
            low_count / len(positions)
            if positions else 0
        ),
        "segment_count": len(runs),
        "has_segment": int(len(runs) > 0),
    })

    s = species_stats[accession]
    s["proteins"] += 1
    s["eligible_positions"] += len(positions)
    s["low_positions"] += low_count
    s["segments"] += len(runs)

    if runs:
        s["proteins_with_segment"] += 1


# --------------------------------------------------
# 7. Write segments
# --------------------------------------------------

with OUT_SEGMENTS.open("w", newline="") as f:

    fields = [
        "orthogroup",
        "accession",
        "protein_id",
        "sequence_id",
        "segment_index",
        "segment_start",
        "segment_end",
        "segment_center",
        "segment_length_codons",
        "codons",
        "weights",
        "species_threshold",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(segment_rows)


# --------------------------------------------------
# 8. Write protein summary
# --------------------------------------------------

with OUT_PROTEINS.open("w", newline="") as f:

    fields = [
        "orthogroup",
        "accession",
        "protein_id",
        "sequence_id",
        "protein_length",
        "eligible_positions",
        "low_adaptation_positions",
        "low_adaptation_fraction",
        "segment_count",
        "has_segment",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(protein_rows)


# --------------------------------------------------
# 9. Species summary
# --------------------------------------------------

species_rows = []

for accession in sorted(species_stats):

    s = species_stats[accession]

    species_rows.append({
        "accession": accession,
        "proteins": s["proteins"],
        "proteins_with_segment": s["proteins_with_segment"],
        "fraction_proteins_with_segment":
            s["proteins_with_segment"] / s["proteins"],
        "eligible_positions": s["eligible_positions"],
        "low_adaptation_positions": s["low_positions"],
        "low_adaptation_fraction":
            s["low_positions"] / s["eligible_positions"],
        "segments": s["segments"],
    })

with OUT_SPECIES.open("w", newline="") as f:

    fields = [
        "accession",
        "proteins",
        "proteins_with_segment",
        "fraction_proteins_with_segment",
        "eligible_positions",
        "low_adaptation_positions",
        "low_adaptation_fraction",
        "segments",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(species_rows)


# --------------------------------------------------
# 10. Summary
# --------------------------------------------------

proteins_with_segment = sum(
    r["has_segment"] for r in protein_rows
)

families_with_segment = len({
    r["orthogroup"]
    for r in segment_rows
})

length_counts = Counter(
    r["segment_length_codons"]
    for r in segment_rows
)

print()
print("Proteins analyzed:", len(protein_rows))
print("Low-adaptation segments:", len(segment_rows))
print(
    "Proteins with >=1 segment:",
    proteins_with_segment
)
print(
    "Families with >=1 segment:",
    families_with_segment
)

print()
print("Segment length distribution:")
for length, count in sorted(length_counts.items()):
    print(length, count)

print()
print("Outputs:")
print(OUT_SEGMENTS)
print(OUT_PROTEINS)
print(OUT_SPECIES)
