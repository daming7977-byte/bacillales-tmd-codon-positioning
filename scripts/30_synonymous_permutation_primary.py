from pathlib import Path
from collections import defaultdict
import argparse
import csv
import random
import re
import statistics
import math


parser = argparse.ArgumentParser()

parser.add_argument(
    "--permutations",
    type=int,
    default=100
)

parser.add_argument(
    "--seed",
    type=int,
    default=20260825
)

args = parser.parse_args()

N_PERM = args.permutations
SEED = args.seed


# --------------------------------------------------
# paths
# --------------------------------------------------

CDS_BASE = Path("data/cds")

WEIGHT_DIR = Path(
    "work/codon/species_weights"
)

THRESHOLD_FILE = Path(
    "work/codon/frozen_thresholds/"
    "species_low_adaptation_thresholds.tsv"
)

MODAL_FILE = Path(
    "work/topology/merged/"
    "modal_tmd_primary_set_diagnostic.tsv"
)

ALL_CLUSTER_MEMBERS = Path(
    "work/topology/merged/"
    "tmd_overlap_cluster_members.tsv"
)

QUALIFIED_CLUSTERS = Path(
    "work/topology/merged/"
    "topology_qualified_tmd_clusters.tsv"
)

OBSERVED_UNITS = Path(
    "work/codon/primary_tmd_assignment/"
    "primary_homologous_tmd_unit_summary.tsv"
)

OUTDIR = Path(
    "work/codon/permutation"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTFILE = OUTDIR / (
    f"primary_synonymous_null_{N_PERM}.tsv"
)


MIN_RUN = 3


# --------------------------------------------------
# genetic code
# --------------------------------------------------

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


# --------------------------------------------------
# FASTA
# --------------------------------------------------

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
# frozen thresholds
# --------------------------------------------------

thresholds = {}

with THRESHOLD_FILE.open() as f:

    r = csv.DictReader(
        f,
        delimiter="\t"
    )

    for x in r:

        thresholds[x["accession"]] = float(
            x["threshold"]
        )


# --------------------------------------------------
# codon weights
# --------------------------------------------------

weights_by_species = {}

for accession in thresholds:

    weights = {}

    path = (
        WEIGHT_DIR
        / f"{accession}_codon_weights.tsv"
    )

    with path.open() as f:

        r = csv.DictReader(
            f,
            delimiter="\t"
        )

        for x in r:

            weights[x["codon"]] = float(
                x["relative_weight"]
            )

    weights_by_species[accession] = weights


# --------------------------------------------------
# modal proteins
# --------------------------------------------------

modal_ok = {}

with MODAL_FILE.open() as f:

    r = csv.DictReader(
        f,
        delimiter="\t"
    )

    for x in r:

        if x["matches_modal"] != "1":
            continue

        modal_ok[x["sequence_id"]] = {
            "orthogroup": x["orthogroup"],
            "accession": x["accession"],
            "protein_id": x["protein_id"],
        }

print(
    "Modal-matching proteins:",
    len(modal_ok)
)


# --------------------------------------------------
# qualified clusters
# --------------------------------------------------

qualified_keys = set()

with QUALIFIED_CLUSTERS.open() as f:

    r = csv.DictReader(
        f,
        delimiter="\t"
    )

    for x in r:

        qualified_keys.add(
            (
                x["orthogroup"],
                int(x["cluster_index"])
            )
        )

print(
    "Frozen qualified clusters:",
    len(qualified_keys)
)


# --------------------------------------------------
# all TMDs of modal proteins
#
# nearest anchor must be found BEFORE filtering
# cluster qualification
# --------------------------------------------------

tmds_by_sequence = defaultdict(list)

cluster_species = defaultdict(set)

with ALL_CLUSTER_MEMBERS.open() as f:

    r = csv.DictReader(
        f,
        delimiter="\t"
    )

    for x in r:

        sid = x["sequence_id"]

        if sid not in modal_ok:
            continue

        start = int(
            x["native_start"]
        )

        end = int(
            x["native_end"]
        )

        rec = {
            "orthogroup":
                x["orthogroup"],

            "cluster_index":
                int(x["cluster_index"]),

            "original_tmd_index":
                int(x["original_tmd_index"]),

            "native_start":
                start,

            "native_end":
                end,

            "native_center":
                (start + end) / 2,
        }

        tmds_by_sequence[sid].append(rec)

        cluster_species[
            (
                x["orthogroup"],
                int(x["cluster_index"])
            )
        ].add(
            x["accession"]
        )


# --------------------------------------------------
# eligible primary clusters:
# frozen QC + >=8 species after modal filter
# --------------------------------------------------

eligible_clusters = {
    key
    for key in qualified_keys
    if len(
        cluster_species.get(
            key,
            set()
        )
    ) >= 8
}

print(
    "Primary clusters with >=8 modal species:",
    len(eligible_clusters)
)


# --------------------------------------------------
# load CDS for 1372 modal proteins
# --------------------------------------------------

cds_indexes = {}

for accession in thresholds:

    files = list(
        (CDS_BASE / accession).rglob(
            "cds_from_genomic.fna"
        )
    )

    if len(files) != 1:
        raise RuntimeError(
            f"{accession}: "
            f"CDS file count = {len(files)}"
        )

    index = {}

    for header, seq in read_fasta(
        files[0]
    ).items():

        m = re.search(
            r"\[protein_id=([^\]]+)\]",
            header
        )

        if m:
            index[m.group(1)] = seq.upper()

    cds_indexes[accession] = index


# --------------------------------------------------
# prepare protein codon records
#
# initiation and terminal stop excluded
# --------------------------------------------------

proteins = {}

for sid, info in modal_ok.items():

    accession = info["accession"]
    protein_id = info["protein_id"]

    cds = cds_indexes[
        accession
    ].get(
        protein_id
    )

    if cds is None:
        raise RuntimeError(
            f"Missing CDS: {sid}"
        )

    codons = [
        cds[i:i+3]
        for i in range(
            0,
            len(cds),
            3
        )
    ]

    if (
        codons
        and codons[-1]
        in {"TAA", "TAG", "TGA"}
    ):
        codons = codons[:-1]

    # initiation codon excluded
    body_codons = codons[1:]

    aa = []

    for c in body_codons:

        a = CODON_TABLE.get(c)

        if a is None or a == "*":
            raise RuntimeError(
                f"Unexpected codon in {sid}: {c}"
            )

        aa.append(a)

    # group codon indices by amino acid
    aa_positions = defaultdict(list)

    for i, a in enumerate(aa):
        aa_positions[a].append(i)

    proteins[sid] = {
        "orthogroup":
            info["orthogroup"],

        "accession":
            accession,

        "protein_id":
            protein_id,

        "codons":
            body_codons,

        "aa_positions":
            dict(aa_positions),

        "tmds":
            tmds_by_sequence[sid],
    }


print(
    "Protein CDS records prepared:",
    len(proteins)
)


# --------------------------------------------------
# observed statistics
# --------------------------------------------------

observed = {}

with OBSERVED_UNITS.open() as f:

    r = csv.DictReader(
        f,
        delimiter="\t"
    )

    for anchor in [
        "start",
        "end",
        "center"
    ]:

        vals = []

        f.seek(0)
        r = csv.DictReader(
            f,
            delimiter="\t"
        )

        for x in r:

            if (
                x["anchor_type"] == anchor
                and x["qualifying_unit"] == "1"
            ):
                vals.append(
                    float(
                        x["sample_variance"]
                    )
                )

        observed[anchor] = {
            "n_units": len(vals),
            "median_variance":
                statistics.median(vals)
                if vals else math.nan,
        }


print()
print("Observed statistics")

for a in ["start", "end", "center"]:

    print(
        a,
        "units=",
        observed[a]["n_units"],
        "median_variance=",
        observed[a]["median_variance"]
    )


# --------------------------------------------------
# helpers
# --------------------------------------------------

anchor_field = {
    "start": "native_start",
    "end": "native_end",
    "center": "native_center",
}


def call_segments(
    codons,
    accession
):

    weights = weights_by_species[
        accession
    ]

    threshold = thresholds[
        accession
    ]

    low = [
        weights[c] <= threshold
        for c in codons
    ]

    segments = []

    run_start = None

    for i, is_low in enumerate(low):

        # body codon index i corresponds
        # biological codon position i + 2
        biological_pos = i + 2

        if is_low:

            if run_start is None:
                run_start = biological_pos

        else:

            if run_start is not None:

                end = biological_pos - 1

                if (
                    end - run_start + 1
                    >= MIN_RUN
                ):
                    segments.append(
                        (
                            run_start,
                            end,
                            (
                                run_start + end
                            ) / 2
                        )
                    )

                run_start = None

    if run_start is not None:

        end = len(codons) + 1

        if (
            end - run_start + 1
            >= MIN_RUN
        ):
            segments.append(
                (
                    run_start,
                    end,
                    (
                        run_start + end
                    ) / 2
                )
            )

    return segments


def permute_synonymous(
    original,
    aa_positions,
    rng
):

    shuffled = list(original)

    for positions in aa_positions.values():

        if len(positions) <= 1:
            continue

        vals = [
            original[i]
            for i in positions
        ]

        rng.shuffle(vals)

        for i, codon in zip(
            positions,
            vals
        ):
            shuffled[i] = codon

    return shuffled


def collect_assignments(
    segments,
    tmds,
    anchor
):

    field = anchor_field[anchor]

    # same protein/cluster:
    # keep closest assigned segment
    retained = {}

    for seg_start, seg_end, center in segments:

        candidates = []

        for tmd in tmds:

            apos = tmd[field]

            rel = (
                center - apos
            )

            candidates.append(
                (
                    abs(rel),
                    tmd[
                        "original_tmd_index"
                    ],
                    tmd[
                        "cluster_index"
                    ],
                    rel,
                    tmd
                )
            )

        if not candidates:
            continue

        candidates.sort(
            key=lambda z: (
                z[0],
                z[1],
                z[2]
            )
        )

        (
            absdist,
            original_tmd_index,
            cluster_index,
            rel,
            tmd
        ) = candidates[0]

        key = (
            tmd["orthogroup"],
            cluster_index
        )

        # IMPORTANT:
        # if nearest TMD cluster is not
        # eligible, do not reassign
        if key not in eligible_clusters:
            continue

        old = retained.get(key)

        candidate = (
            absdist,
            center,
            rel
        )

        if (
            old is None
            or candidate[:2] < old[:2]
        ):
            retained[key] = candidate

    return {
        key: x[2]
        for key, x in retained.items()
    }


# --------------------------------------------------
# permutation
# --------------------------------------------------

rng = random.Random(SEED)

output_rows = []

for perm in range(
    1,
    N_PERM + 1
):

    unit_positions = {
        "start":
            defaultdict(list),

        "end":
            defaultdict(list),

        "center":
            defaultdict(list),
    }

    # ----------------------------------------------
    # shuffle each protein independently
    # ----------------------------------------------

    for sid, p in proteins.items():

        shuffled = permute_synonymous(
            p["codons"],
            p["aa_positions"],
            rng
        )

        segments = call_segments(
            shuffled,
            p["accession"]
        )

        if not segments:
            continue

        for anchor in [
            "start",
            "end",
            "center"
        ]:

            assigned = collect_assignments(
                segments,
                p["tmds"],
                anchor
            )

            for key, rel in assigned.items():

                # one protein per species per OG,
                # so one value = one species
                unit_positions[
                    anchor
                ][key].append(rel)


    # ----------------------------------------------
    # summarize replicate
    # ----------------------------------------------

    row = {
        "permutation":
            perm
    }

    for anchor in [
        "start",
        "end",
        "center"
    ]:

        variances = []

        for key, vals in (
            unit_positions[
                anchor
            ].items()
        ):

            if len(vals) < 3:
                continue

            variances.append(
                statistics.variance(vals)
            )

        row[
            f"{anchor}_qualifying_units"
        ] = len(variances)

        row[
            f"{anchor}_median_variance"
        ] = (
            statistics.median(
                variances
            )
            if variances
            else ""
        )

    output_rows.append(row)

    if (
        perm <= 10
        or perm % 25 == 0
        or perm == N_PERM
    ):

        print(
            f"Permutation {perm}/{N_PERM}",
            "start_units=",
            row[
                "start_qualifying_units"
            ],
            "start_median=",
            row[
                "start_median_variance"
            ],
            "end_units=",
            row[
                "end_qualifying_units"
            ],
            "end_median=",
            row[
                "end_median_variance"
            ]
        )


# --------------------------------------------------
# write
# --------------------------------------------------

fields = [
    "permutation",

    "start_qualifying_units",
    "start_median_variance",

    "end_qualifying_units",
    "end_median_variance",

    "center_qualifying_units",
    "center_median_variance",
]

with OUTFILE.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(output_rows)


# --------------------------------------------------
# final comparison
# --------------------------------------------------

print()
print("NULL SUMMARY")

for anchor in [
    "start",
    "end",
    "center"
]:

    vals = [
        float(
            x[
                f"{anchor}_median_variance"
            ]
        )
        for x in output_rows
        if x[
            f"{anchor}_median_variance"
        ] != ""
    ]

    units = [
        int(
            x[
                f"{anchor}_qualifying_units"
            ]
        )
        for x in output_rows
    ]

    obs = observed[
        anchor
    ]["median_variance"]

    k = sum(
        x <= obs
        for x in vals
    )

    p = (
        (k + 1)
        /
        (len(vals) + 1)
        if vals
        else math.nan
    )

    print()
    print("Anchor:", anchor)

    print(
        "  observed units:",
        observed[
            anchor
        ]["n_units"]
    )

    print(
        "  observed median variance:",
        obs
    )

    print(
        "  null qualifying units median:",
        statistics.median(units)
    )

    print(
        "  null qualifying units range:",
        min(units),
        max(units)
    )

    print(
        "  null median-of-medians:",
        statistics.median(vals)
    )

    print(
        "  null range:",
        min(vals),
        max(vals)
    )

    print(
        "  null <= observed:",
        k,
        "/",
        len(vals)
    )

    print(
        "  empirical one-sided P:",
        p
    )


print()
print("Output:", OUTFILE)
