from pathlib import Path
from collections import defaultdict

PANEL = Path("config/panel_conservative_10.tsv")
GROUPS = Path(
    "work/orthology/strict_groups_conservative10/"
    "strict_one_to_one_groups.tsv"
)

OUTDIR = Path("work/topology")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_FASTA = OUTDIR / "bacillales_ge8of10_proteins.faa"
OUT_GROUPS = OUTDIR / "bacillales_ge8of10_groups.tsv"

# --------------------------------------------------
# Read panel accessions
# --------------------------------------------------
accessions = []
with PANEL.open() as f:
    header = next(f)
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        accessions.append(fields[0])

print("Panel accessions:", len(accessions))

# --------------------------------------------------
# Load all proteomes
# key = accession|protein_id
# --------------------------------------------------
sequences = {}

for acc in accessions:
    faa = Path(
        f"data/refseq/{acc}/ncbi_dataset/data/{acc}/protein.faa"
    )

    if not faa.exists():
        raise FileNotFoundError(faa)

    current_id = None
    seq_chunks = []

    with faa.open() as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if current_id is not None:
                    sequences[f"{acc}|{current_id}"] = "".join(seq_chunks)

                current_id = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

        if current_id is not None:
            sequences[f"{acc}|{current_id}"] = "".join(seq_chunks)

print("Proteins loaded:", len(sequences))

# --------------------------------------------------
# Read strict orthogroups and retain >=8/10
# --------------------------------------------------
selected_groups = []

with GROUPS.open() as f:
    header = next(f).rstrip("\n").split("\t")

    for line in f:
        fields = line.rstrip("\n").split("\t")

        orthogroup = fields[0]
        species_count = int(fields[1])
        members = fields[2].split(";")

        if species_count >= 8:
            selected_groups.append(
                (orthogroup, species_count, members)
            )

print("Orthogroups >=8/10:", len(selected_groups))

# --------------------------------------------------
# Validate and write
# --------------------------------------------------
missing = []
protein_count = 0

with OUT_FASTA.open("w") as fout, OUT_GROUPS.open("w") as gout:

    gout.write("orthogroup\tspecies_count\tmember_count\tmembers\n")

    for orthogroup, species_count, members in selected_groups:

        gout.write(
            f"{orthogroup}\t{species_count}\t{len(members)}\t"
            + ",".join(members)
            + "\n"
        )

        for member in members:
            if member not in sequences:
                missing.append((orthogroup, member))
                continue

            seq = sequences[member]

            fout.write(
                f">{orthogroup}|{member}\n{seq}\n"
            )

            protein_count += 1

print("Proteins written:", protein_count)
print("Missing proteins:", len(missing))
print("FASTA:", OUT_FASTA)
print("Group table:", OUT_GROUPS)

if missing:
    print("\nFirst missing entries:")
    for x in missing[:20]:
        print(x)
