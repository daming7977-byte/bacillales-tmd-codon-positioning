from pathlib import Path
from collections import Counter
import csv
import re
import math

CDS_BASE = Path("data/cds")
WEIGHT_DIR = Path("work/codon/species_weights")

ACCESSIONS = [
    "GCF_000009045.1",
    "GCA_000015785.2",
    "GCA_000025805.1",
    "GCF_000011245.1",
    "GCF_034478925.1",
    "GCF_006094295.1",
    "GCF_000017885.4",
    "GCF_000009785.1",
    "GCF_034479305.1",
    "GCF_005671335.1",
]

def read_fasta(path):
    records = []
    header = None
    seq = []

    with path.open() as f:
        for line in f:
            line = line.rstrip()

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())

        if header is not None:
            records.append((header, "".join(seq)))

    return records

def lower_decile(values):
    values = sorted(values)

    if not values:
        return None

    # nearest-rank style empirical 10th percentile:
    # choose value at ceil(0.10*N), 1-based
    rank = max(1, math.ceil(0.10 * len(values)))
    return values[rank - 1]

print(
    "accession\t"
    "unique61_threshold\t"
    "position_weighted_threshold\t"
    "unique61_n\t"
    "position_n"
)

for acc in ACCESSIONS:

    # -----------------------------
    # load codon weights
    # -----------------------------
    weights = {}

    wf = WEIGHT_DIR / f"{acc}_codon_weights.tsv"

    with wf.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            weights[row["codon"]] = float(row["relative_weight"])

    # -----------------------------
    # A: threshold across unique sense codons
    # -----------------------------
    unique_values = list(weights.values())
    unique_thr = lower_decile(unique_values)

    # -----------------------------
    # B: threshold across actual codon positions
    #    excluding start and terminal stop
    # -----------------------------
    cds_files = list(
        (CDS_BASE / acc).rglob("cds_from_genomic.fna")
    )

    if len(cds_files) != 1:
        raise RuntimeError(
            f"{acc}: expected one CDS file, found {len(cds_files)}"
        )

    position_values = []

    for header, seq in read_fasta(cds_files[0]):

        seq = seq.upper().replace("U", "T")

        if len(seq) % 3 != 0:
            continue

        if any(base not in "ACGT" for base in seq):
            continue

        codons = [
            seq[i:i+3]
            for i in range(0, len(seq), 3)
        ]

        # exclude initiation codon
        codons = codons[1:]

        # exclude terminal stop if present
        if codons and codons[-1] in {"TAA", "TAG", "TGA"}:
            codons = codons[:-1]

        for codon in codons:
            if codon in weights:
                position_values.append(weights[codon])

    pos_thr = lower_decile(position_values)

    print(
        f"{acc}\t"
        f"{unique_thr:.6f}\t"
        f"{pos_thr:.6f}\t"
        f"{len(unique_values)}\t"
        f"{len(position_values)}"
    )
