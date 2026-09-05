from pathlib import Path
from collections import Counter, defaultdict
import csv
import re

CDS_BASE = Path("data/cds")
OUTDIR = Path("work/codon/species_weights")
OUTDIR.mkdir(parents=True, exist_ok=True)

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


summary_rows = []

for accession in ACCESSIONS:

    matches = list(
        (CDS_BASE / accession).rglob("cds_from_genomic.fna")
    )

    if len(matches) != 1:
        raise RuntimeError(
            f"{accession}: expected 1 CDS file, found {len(matches)}"
        )

    cds_file = matches[0]

    records = read_fasta(cds_file)

    codon_counts = Counter()

    total_records = 0
    usable_records = 0
    excluded_bad_length = 0
    excluded_ambiguous = 0
    excluded_internal_stop = 0

    for header, seq in records:

        total_records += 1

        seq = seq.upper().replace("U", "T")

        if len(seq) % 3 != 0:
            excluded_bad_length += 1
            continue

        if any(base not in "ACGT" for base in seq):
            excluded_ambiguous += 1
            continue

        codons = [
            seq[i:i+3]
            for i in range(0, len(seq), 3)
        ]

        aas = [
            CODON_TABLE.get(codon, "X")
            for codon in codons
        ]

        # terminal stop is allowed
        internal = aas[:-1] if aas and aas[-1] == "*" else aas

        if "*" in internal:
            excluded_internal_stop += 1
            continue

        usable_records += 1

        # Stop codons excluded from codon usage
        for codon, aa in zip(codons, aas):
            if aa == "*":
                continue

            codon_counts[codon] += 1


    # --------------------------------------------------
    # synonymous codon groups
    # --------------------------------------------------

    aa_to_codons = defaultdict(list)

    for codon, aa in CODON_TABLE.items():
        if aa == "*":
            continue

        aa_to_codons[aa].append(codon)


    weight_rows = []

    for aa in sorted(aa_to_codons):

        synonymous = aa_to_codons[aa]

        max_count = max(
            codon_counts[c] for c in synonymous
        )

        for codon in sorted(synonymous):

            count = codon_counts[codon]

            if max_count > 0:
                weight = count / max_count
            else:
                weight = 0.0

            weight_rows.append({
                "accession": accession,
                "amino_acid": aa,
                "codon": codon,
                "count": count,
                "relative_weight": weight,
            })


    outfile = OUTDIR / f"{accession}_codon_weights.tsv"

    with outfile.open("w", newline="") as f:

        fields = [
            "accession",
            "amino_acid",
            "codon",
            "count",
            "relative_weight",
        ]

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(weight_rows)


    summary_rows.append({
        "accession": accession,
        "total_cds_records": total_records,
        "usable_cds_records": usable_records,
        "excluded_bad_length": excluded_bad_length,
        "excluded_ambiguous": excluded_ambiguous,
        "excluded_internal_stop": excluded_internal_stop,
        "total_nonstop_codons": sum(codon_counts.values()),
    })

    print(
        accession,
        "total=", total_records,
        "usable=", usable_records,
        "codons=", sum(codon_counts.values())
    )


summary_file = OUTDIR / "species_codon_weight_summary.tsv"

with summary_file.open("w", newline="") as f:

    fields = [
        "accession",
        "total_cds_records",
        "usable_cds_records",
        "excluded_bad_length",
        "excluded_ambiguous",
        "excluded_internal_stop",
        "total_nonstop_codons",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(summary_rows)


print()
print("Species processed:", len(summary_rows))
print("Output directory:", OUTDIR)
print("Summary:", summary_file)
