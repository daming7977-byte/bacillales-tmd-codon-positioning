from pathlib import Path
from collections import defaultdict
import csv

SEGMENTS = Path(
    "work/topology/merged/deeptmhmm_all_segments.tsv"
)

ALNDIR = Path(
    "work/topology/candidate_alignments"
)

OUTFILE = Path(
    "work/topology/merged/tmd_alignment_coordinates.tsv"
)


def read_fasta(path):
    records = {}
    seq_id = None
    seq = []

    with path.open() as f:
        for line in f:
            line = line.rstrip()

            if line.startswith(">"):
                if seq_id is not None:
                    records[seq_id] = "".join(seq)

                seq_id = line[1:].split()[0]
                seq = []

            else:
                seq.append(line)

        if seq_id is not None:
            records[seq_id] = "".join(seq)

    return records


# --------------------------------------------------
# 1. Load DeepTMHMM TMhelix coordinates
# --------------------------------------------------

tmds = defaultdict(list)

with SEGMENTS.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        if row["segment_label"] != "TMhelix":
            continue

        seq_id = (
            f'{row["orthogroup"]}|'
            f'{row["accession"]}|'
            f'{row["protein_id"]}'
        )

        start = int(row["start"])
        end = int(row["end"])

        tmds[seq_id].append((start, end))

# Ensure N -> C order
for seq_id in tmds:
    tmds[seq_id].sort()

print("Proteins with >=1 TMhelix:", len(tmds))
print(
    "TMhelix segments loaded:",
    sum(len(x) for x in tmds.values())
)

# --------------------------------------------------
# 2. Map native residue coordinates to alignment cols
# --------------------------------------------------

rows = []

alignment_files = sorted(
    ALNDIR.glob("orthogroup_*.aln.faa"),
    key=lambda p: int(p.stem.split("_")[1].split(".")[0])
)

mapping_errors = []
sequence_count = 0
mapped_tmd_count = 0

for alnfile in alignment_files:

    orthogroup = alnfile.name.replace(".aln.faa", "")
    aligned = read_fasta(alnfile)

    for seq_id, alnseq in aligned.items():

        sequence_count += 1

        parts = seq_id.split("|", 2)

        if len(parts) != 3:
            mapping_errors.append(
                (seq_id, "unexpected sequence ID")
            )
            continue

        og, accession, protein_id = parts

        # residue number (1-based) -> alignment column (1-based)
        residue_to_col = {}

        residue_number = 0

        for col, aa in enumerate(alnseq, start=1):
            if aa != "-":
                residue_number += 1
                residue_to_col[residue_number] = col

        seq_tmds = tmds.get(seq_id, [])

        for tmd_index, (native_start, native_end) in enumerate(
            seq_tmds,
            start=1
        ):

            if (
                native_start not in residue_to_col
                or native_end not in residue_to_col
            ):
                mapping_errors.append(
                    (
                        seq_id,
                        f"TMD {tmd_index}: "
                        f"{native_start}-{native_end} "
                        f"outside sequence length {residue_number}"
                    )
                )
                continue

            aln_start = residue_to_col[native_start]
            aln_end = residue_to_col[native_end]

            rows.append({
                "orthogroup": og,
                "accession": accession,
                "protein_id": protein_id,
                "sequence_id": seq_id,
                "tmd_index": tmd_index,
                "native_start": native_start,
                "native_end": native_end,
                "native_length": native_end - native_start + 1,
                "alignment_start": aln_start,
                "alignment_end": aln_end,
                "alignment_span": aln_end - aln_start + 1,
            })

            mapped_tmd_count += 1

# --------------------------------------------------
# 3. Output
# --------------------------------------------------

with OUTFILE.open("w", newline="") as f:

    fieldnames = [
        "orthogroup",
        "accession",
        "protein_id",
        "sequence_id",
        "tmd_index",
        "native_start",
        "native_end",
        "native_length",
        "alignment_start",
        "alignment_end",
        "alignment_span",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

print("Alignment files:", len(alignment_files))
print("Aligned sequences examined:", sequence_count)
print("Mapped TMhelix segments:", mapped_tmd_count)
print("Mapping errors:", len(mapping_errors))
print("Output:", OUTFILE)

if mapping_errors:
    print("\nFirst mapping errors:")
    for x in mapping_errors[:20]:
        print(x)
