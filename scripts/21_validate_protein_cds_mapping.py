from pathlib import Path
from collections import defaultdict
import csv
import re

PROTEIN_FASTA = Path(
    "work/topology/bacillales_ge8of10_proteins.faa"
)

FAMILY_TABLE = Path(
    "work/topology/merged/topology_qualified_families.tsv"
)

CDS_BASE = Path("data/cds")

OUT = Path(
    "work/codon/protein_cds_validation.tsv"
)
OUT.parent.mkdir(parents=True, exist_ok=True)


CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


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


def translate(seq):
    seq = seq.upper().replace("U", "T")

    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]

        if len(codon) != 3:
            break

        if any(x not in "ACGT" for x in codon):
            aa.append("X")
        else:
            aa.append(CODON_TABLE.get(codon, "X"))

    return "".join(aa)


# --------------------------------------------------
# 1. Qualified orthogroups
# --------------------------------------------------

qualified = set()

with FAMILY_TABLE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        if row["topology_qualified"] == "1":
            qualified.add(row["orthogroup"])

print("Qualified families:", len(qualified))


# --------------------------------------------------
# 2. Load target proteins
# --------------------------------------------------

protein_records = read_fasta(PROTEIN_FASTA)

targets = {}

for header, seq in protein_records.items():

    seq_id = header.split()[0]
    parts = seq_id.split("|", 2)

    if len(parts) != 3:
        continue

    og, accession, protein_id = parts

    if og in qualified:
        targets[seq_id] = {
            "orthogroup": og,
            "accession": accession,
            "protein_id": protein_id,
            "protein_seq": seq,
        }

print("Target proteins:", len(targets))


# --------------------------------------------------
# 3. Build protein_id -> CDS mapping
# --------------------------------------------------

cds_index = defaultdict(dict)

for accession in sorted(set(x["accession"] for x in targets.values())):

    matches = list(
        (CDS_BASE / accession).rglob("cds_from_genomic.fna")
    )

    if len(matches) != 1:
        raise RuntimeError(
            f"{accession}: expected exactly 1 CDS file, found {len(matches)}"
        )

    cds_file = matches[0]

    print("Loading CDS:", accession, cds_file)

    records = read_fasta(cds_file)

    for header, seq in records.items():

        # NCBI CDS FASTA headers commonly contain:
        # [protein_id=WP_....]
        m = re.search(r"\[protein_id=([^\]]+)\]", header)

        if not m:
            continue

        protein_id = m.group(1)

        cds_index[accession][protein_id] = seq


# --------------------------------------------------
# 4. Validate
# --------------------------------------------------

rows = []

for seq_id, rec in targets.items():

    og = rec["orthogroup"]
    accession = rec["accession"]
    protein_id = rec["protein_id"]
    protein = rec["protein_seq"]

    cds = cds_index.get(accession, {}).get(protein_id)

    if cds is None:
        rows.append({
            "orthogroup": og,
            "accession": accession,
            "protein_id": protein_id,
            "protein_length": len(protein),
            "cds_length": "",
            "translation_length": "",
            "status": "CDS_NOT_FOUND",
        })
        continue

    translated = translate(cds)

    # remove terminal stop if present
    translated_no_stop = (
        translated[:-1]
        if translated.endswith("*")
        else translated
    )

    if translated_no_stop == protein:
        status = "EXACT_MATCH"

    elif (
        len(translated_no_stop) == len(protein)
        and len(protein) > 0
        and protein[0] == "M"
        and translated_no_stop[1:] == protein[1:]
    ):
        status = "VALID_ALTERNATIVE_START"

    elif len(translated_no_stop) == len(protein):
        status = "SAME_LENGTH_MISMATCH"

    else:
        status = "LENGTH_MISMATCH"

    rows.append({
        "orthogroup": og,
        "accession": accession,
        "protein_id": protein_id,
        "protein_length": len(protein),
        "cds_length": len(cds),
        "translation_length": len(translated_no_stop),
        "status": status,
    })


# --------------------------------------------------
# 5. Output
# --------------------------------------------------

with OUT.open("w", newline="") as f:
    fieldnames = [
        "orthogroup",
        "accession",
        "protein_id",
        "protein_length",
        "cds_length",
        "translation_length",
        "status",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(rows)


# --------------------------------------------------
# 6. Summary
# --------------------------------------------------

from collections import Counter

counts = Counter(r["status"] for r in rows)

print()
print("Validation summary:")
for k, v in sorted(counts.items()):
    print(f"{k}\t{v}")

print()
print("Total validated targets:", len(rows))
print("Output:", OUT)
