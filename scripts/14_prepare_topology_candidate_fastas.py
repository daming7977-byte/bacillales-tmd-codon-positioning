from pathlib import Path
from collections import defaultdict
import csv

SUMMARY = Path(
    "work/topology/merged/orthogroup_topology_summary.tsv"
)

FASTA = Path(
    "work/topology/bacillales_ge8of10_proteins.faa"
)

OUTDIR = Path(
    "work/topology/candidate_fastas"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# 1. Select candidate orthogroups
#    Criterion: modal nTMD >= 3
# --------------------------------------------------

candidate_ogs = set()

with SUMMARY.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        if int(float(row["modal_n_tmd"])) >= 3:
            candidate_ogs.add(row["orthogroup"])

print("Candidate orthogroups:", len(candidate_ogs))

# --------------------------------------------------
# 2. Read master FASTA
# --------------------------------------------------

seqs = {}

current_id = None
current_seq = []

with FASTA.open() as f:
    for line in f:
        line = line.rstrip()

        if line.startswith(">"):
            if current_id is not None:
                seqs[current_id] = "".join(current_seq)

            current_id = line[1:].split()[0]
            current_seq = []

        else:
            current_seq.append(line)

    if current_id is not None:
        seqs[current_id] = "".join(current_seq)

print("Sequences loaded:", len(seqs))

# --------------------------------------------------
# 3. Group candidate sequences
# --------------------------------------------------

groups = defaultdict(list)

for seq_id, seq in seqs.items():

    parts = seq_id.split("|", 2)

    if len(parts) != 3:
        raise ValueError(f"Unexpected FASTA ID: {seq_id}")

    orthogroup = parts[0]

    if orthogroup in candidate_ogs:
        groups[orthogroup].append((seq_id, seq))

# --------------------------------------------------
# 4. Write one FASTA per orthogroup
# --------------------------------------------------

written_sequences = 0

for og in sorted(
    groups,
    key=lambda x: int(x.split("_")[-1])
):

    outfile = OUTDIR / f"{og}.faa"

    with outfile.open("w") as out:
        for seq_id, seq in groups[og]:
            out.write(f">{seq_id}\n")

            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

            written_sequences += 1

print("Candidate FASTA files written:", len(groups))
print("Candidate sequences written:", written_sequences)
print("Output directory:", OUTDIR)

missing = candidate_ogs - set(groups)

print("Candidate orthogroups missing from FASTA:", len(missing))

if missing:
    for og in sorted(missing):
        print("MISSING:", og)
