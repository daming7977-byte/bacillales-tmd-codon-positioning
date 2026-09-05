from pathlib import Path
from collections import defaultdict
import csv
import re
import math

CDS_BASE = Path("data/cds")
WEIGHT_DIR = Path("work/codon/species_weights")

FAMILY_TABLE = Path(
    "work/topology/merged/topology_qualified_families.tsv"
)

PROTEIN_FASTA = Path(
    "work/topology/bacillales_ge8of10_proteins.faa"
)

OUTDIR = Path("work/codon/frozen_thresholds")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "species_low_adaptation_thresholds.tsv"


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
# 1. Qualified orthogroups
# --------------------------------------------------

qualified = set()

with FAMILY_TABLE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        if row["topology_qualified"] == "1":
            qualified.add(row["orthogroup"])

print("Qualified membrane families:", len(qualified))


# --------------------------------------------------
# 2. Target protein IDs by species
# --------------------------------------------------

targets = defaultdict(set)

for header in read_fasta(PROTEIN_FASTA):

    seq_id = header.split()[0]
    parts = seq_id.split("|", 2)

    if len(parts) != 3:
        continue

    og, accession, protein_id = parts

    if og in qualified:
        targets[accession].add(protein_id)

print(
    "Target membrane proteins:",
    sum(len(x) for x in targets.values())
)


# --------------------------------------------------
# 3. Process each species
# --------------------------------------------------

rows = []

for accession in sorted(targets):

    # Load codon weights
    weights = {}

    wf = WEIGHT_DIR / f"{accession}_codon_weights.tsv"

    with wf.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            weights[row["codon"]] = float(
                row["relative_weight"]
            )

    # Load CDS
    matches = list(
        (CDS_BASE / accession).rglob(
            "cds_from_genomic.fna"
        )
    )

    if len(matches) != 1:
        raise RuntimeError(
            f"{accession}: CDS file count = {len(matches)}"
        )

    cds_records = read_fasta(matches[0])

    protein_to_cds = {}

    for header, seq in cds_records.items():

        m = re.search(
            r"\[protein_id=([^\]]+)\]",
            header
        )

        if m:
            protein_to_cds[m.group(1)] = seq.upper()

    position_weights = []

    proteins_used = 0
    missing = []

    for protein_id in sorted(targets[accession]):

        cds = protein_to_cds.get(protein_id)

        if cds is None:
            missing.append(protein_id)
            continue

        codons = [
            cds[i:i+3]
            for i in range(0, len(cds), 3)
        ]

        # Exclude annotated initiation codon
        codons = codons[1:]

        # Exclude terminal stop
        if codons and codons[-1] in {
            "TAA", "TAG", "TGA"
        }:
            codons = codons[:-1]

        for codon in codons:

            if codon not in weights:
                raise RuntimeError(
                    f"{accession} {protein_id}: "
                    f"unknown codon {codon}"
                )

            position_weights.append(weights[codon])

        proteins_used += 1

    if missing:
        raise RuntimeError(
            f"{accession}: missing {len(missing)} CDS"
        )

    position_weights.sort()

    n = len(position_weights)

    # Empirical bottom decile:
    # observed value at ceil(0.10 * n)
    rank = max(1, math.ceil(0.10 * n))
    threshold = position_weights[rank - 1]

    below = sum(
        x < threshold for x in position_weights
    )

    at = sum(
        x == threshold for x in position_weights
    )

    at_or_below = below + at

    rows.append({
        "accession": accession,
        "proteins_used": proteins_used,
        "eligible_codon_positions": n,
        "bottom_decile_rank": rank,
        "threshold": threshold,
        "positions_below_threshold": below,
        "positions_equal_threshold": at,
        "positions_at_or_below_threshold": at_or_below,
        "fraction_at_or_below_threshold": at_or_below / n,
    })

    print(
        accession,
        "proteins=", proteins_used,
        "positions=", n,
        "threshold=", f"{threshold:.6f}",
        "classified_low=", f"{at_or_below/n:.4f}"
    )


with OUTFILE.open("w", newline="") as f:

    fields = [
        "accession",
        "proteins_used",
        "eligible_codon_positions",
        "bottom_decile_rank",
        "threshold",
        "positions_below_threshold",
        "positions_equal_threshold",
        "positions_at_or_below_threshold",
        "fraction_at_or_below_threshold",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)


print()
print("Species:", len(rows))
print("Frozen thresholds:", OUTFILE)
