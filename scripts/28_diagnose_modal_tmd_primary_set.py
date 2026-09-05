from pathlib import Path
from collections import defaultdict, Counter
import csv

TOPOLOGY_SUMMARY = Path(
    "work/topology/merged/orthogroup_topology_summary.tsv"
)

QUALIFIED_FAMILIES = Path(
    "work/topology/merged/topology_qualified_families.tsv"
)

TMD_COORDS = Path(
    "work/topology/merged/tmd_alignment_coordinates.tsv"
)

PROTEIN_FASTA = Path(
    "work/topology/bacillales_ge8of10_proteins.faa"
)

OUT = Path(
    "work/topology/merged/modal_tmd_primary_set_diagnostic.tsv"
)


def read_fasta_ids(path):
    ids = []

    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                ids.append(
                    line[1:].strip().split()[0]
                )

    return ids


# --------------------------------------------------
# 1. modal TMD count per orthogroup
# --------------------------------------------------

modal = {}

with TOPOLOGY_SUMMARY.open() as f:
    r = csv.DictReader(f, delimiter="\t")

    for x in r:
        modal[x["orthogroup"]] = int(
            x["modal_n_tmd"]
        )


# --------------------------------------------------
# 2. topology-qualified families
# --------------------------------------------------

qualified = set()

with QUALIFIED_FAMILIES.open() as f:
    r = csv.DictReader(f, delimiter="\t")

    for x in r:
        if x["topology_qualified"] == "1":
            qualified.add(x["orthogroup"])

print("Topology-qualified families:", len(qualified))


# --------------------------------------------------
# 3. actual nTMD per sequence
# --------------------------------------------------

ntmd = Counter()

with TMD_COORDS.open() as f:
    r = csv.DictReader(f, delimiter="\t")

    for x in r:
        ntmd[x["sequence_id"]] += 1


# --------------------------------------------------
# 4. qualified-family target proteins
# --------------------------------------------------

rows = []

for sid in read_fasta_ids(PROTEIN_FASTA):

    parts = sid.split("|", 2)

    if len(parts) != 3:
        continue

    og, accession, protein_id = parts

    if og not in qualified:
        continue

    observed = ntmd.get(sid, 0)
    expected = modal[og]

    rows.append({
        "orthogroup": og,
        "accession": accession,
        "protein_id": protein_id,
        "sequence_id": sid,
        "observed_n_tmd": observed,
        "modal_n_tmd": expected,
        "matches_modal": int(
            observed == expected
        ),
    })


# --------------------------------------------------
# 5. family-level retained species
# --------------------------------------------------

family_total = Counter()
family_retained = Counter()

for x in rows:
    og = x["orthogroup"]

    family_total[og] += 1

    if x["matches_modal"]:
        family_retained[og] += 1


with OUT.open("w", newline="") as f:

    fields = [
        "orthogroup",
        "accession",
        "protein_id",
        "sequence_id",
        "observed_n_tmd",
        "modal_n_tmd",
        "matches_modal",
    ]

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(rows)


# --------------------------------------------------
# 6. summaries
# --------------------------------------------------

retained = [
    x for x in rows
    if x["matches_modal"] == 1
]

excluded = [
    x for x in rows
    if x["matches_modal"] == 0
]


print("Qualified-family proteins:", len(rows))
print("Modal-matching proteins:", len(retained))
print("Excluded by modal filter:", len(excluded))

if rows:
    print(
        "Fraction retained:",
        f"{len(retained)/len(rows):.4f}"
    )


print()
print("Families by retained species count:")

dist = Counter(
    family_retained[og]
    for og in qualified
)

for n in sorted(dist):
    print(n, dist[n])


print()
print(
    "Families retaining >=8 modal-matching species:",
    sum(
        family_retained[og] >= 8
        for og in qualified
    )
)

print(
    "Families retaining <8 modal-matching species:",
    sum(
        family_retained[og] < 8
        for og in qualified
    )
)


print()
print("Excluded proteins by nTMD difference:")

diffs = Counter(
    x["observed_n_tmd"] - x["modal_n_tmd"]
    for x in excluded
)

for diff in sorted(diffs):
    print(
        f"{diff:+d}",
        diffs[diff]
    )


print()
print("Families losing most proteins:")

losses = []

for og in qualified:

    total = family_total[og]
    keep = family_retained[og]

    losses.append(
        (total - keep, og, total, keep, modal[og])
    )

for lost, og, total, keep, m in sorted(
    losses,
    reverse=True
)[:20]:

    if lost == 0:
        break

    print(
        og,
        "total=", total,
        "retained=", keep,
        "lost=", lost,
        "modal_nTMD=", m
    )


print()
print("Output:", OUT)
