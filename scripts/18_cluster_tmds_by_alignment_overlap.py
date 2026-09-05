from pathlib import Path
from collections import defaultdict, Counter
import csv
import statistics

INFILE = Path(
    "work/topology/merged/tmd_alignment_coordinates.tsv"
)

FAMILY_SUMMARY = Path(
    "work/topology/merged/orthogroup_topology_summary.tsv"
)

OUT_CLUSTERS = Path(
    "work/topology/merged/tmd_overlap_clusters.tsv"
)

OUT_MEMBERS = Path(
    "work/topology/merged/tmd_overlap_cluster_members.tsv"
)

MIN_OVERLAP_FRACTION = 0.50


# --------------------------------------------------
# Utility
# --------------------------------------------------

def overlap_fraction(a_start, a_end, b_start, b_end):
    """
    Fraction of the shorter interval covered by the overlap.
    Coordinates are inclusive.
    """
    overlap = max(
        0,
        min(a_end, b_end) - max(a_start, b_start) + 1
    )

    len_a = a_end - a_start + 1
    len_b = b_end - b_start + 1

    return overlap / min(len_a, len_b)


def median_abs_deviation(values):
    med = statistics.median(values)
    return statistics.median(
        abs(x - med) for x in values
    )


# --------------------------------------------------
# Load family sizes
# --------------------------------------------------

family_size = {}

with FAMILY_SUMMARY.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        family_size[row["orthogroup"]] = int(row["member_count"])


# --------------------------------------------------
# Load TMDs
# --------------------------------------------------

families = defaultdict(list)

with INFILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for i, row in enumerate(reader):

        row["alignment_start"] = int(row["alignment_start"])
        row["alignment_end"] = int(row["alignment_end"])
        row["native_start"] = int(row["native_start"])
        row["native_end"] = int(row["native_end"])
        row["tmd_index"] = int(row["tmd_index"])

        row["_node_id"] = i

        families[row["orthogroup"]].append(row)


# --------------------------------------------------
# Connected components by interval overlap
# --------------------------------------------------

cluster_rows = []
member_rows = []

for og, tmds in families.items():

    n = len(tmds)

    # Union-find
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra = find(a)
        rb = find(b)

        if ra != rb:
            parent[rb] = ra

    # Compare TMDs from different proteins only
    for i in range(n):
        a = tmds[i]

        for j in range(i + 1, n):
            b = tmds[j]

            if a["sequence_id"] == b["sequence_id"]:
                continue

            frac = overlap_fraction(
                a["alignment_start"],
                a["alignment_end"],
                b["alignment_start"],
                b["alignment_end"],
            )

            if frac >= MIN_OVERLAP_FRACTION:
                union(i, j)

    components = defaultdict(list)

    for i in range(n):
        components[find(i)].append(tmds[i])

    # Sort clusters in N -> C alignment order
    comps = list(components.values())

    comps.sort(
        key=lambda xs: statistics.median(
            x["alignment_start"] for x in xs
        )
    )

    for cluster_index, members in enumerate(comps, start=1):

        seq_counts = Counter(
            x["sequence_id"] for x in members
        )

        unique_sequences = len(seq_counts)

        duplicate_sequences = sum(
            count > 1 for count in seq_counts.values()
        )

        starts = [x["alignment_start"] for x in members]
        ends = [x["alignment_end"] for x in members]

        fam_n = family_size[og]

        cluster_rows.append({
            "orthogroup": og,
            "cluster_index": cluster_index,
            "family_member_count": fam_n,
            "tmd_count": len(members),
            "support_n": unique_sequences,
            "support_fraction": unique_sequences / fam_n,
            "proteins_with_multiple_tmds_in_cluster": duplicate_sequences,
            "median_alignment_start": statistics.median(starts),
            "median_alignment_end": statistics.median(ends),
            "start_min": min(starts),
            "start_max": max(starts),
            "start_range": max(starts) - min(starts),
            "start_mad": median_abs_deviation(starts),
            "end_min": min(ends),
            "end_max": max(ends),
            "end_range": max(ends) - min(ends),
            "end_mad": median_abs_deviation(ends),
        })

        for x in members:
            member_rows.append({
                "orthogroup": og,
                "cluster_index": cluster_index,
                "sequence_id": x["sequence_id"],
                "accession": x["accession"],
                "protein_id": x["protein_id"],
                "original_tmd_index": x["tmd_index"],
                "alignment_start": x["alignment_start"],
                "alignment_end": x["alignment_end"],
                "native_start": x["native_start"],
                "native_end": x["native_end"],
            })


# --------------------------------------------------
# Write outputs
# --------------------------------------------------

with OUT_CLUSTERS.open("w", newline="") as f:

    fields = [
        "orthogroup",
        "cluster_index",
        "family_member_count",
        "tmd_count",
        "support_n",
        "support_fraction",
        "proteins_with_multiple_tmds_in_cluster",
        "median_alignment_start",
        "median_alignment_end",
        "start_min",
        "start_max",
        "start_range",
        "start_mad",
        "end_min",
        "end_max",
        "end_range",
        "end_mad",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(cluster_rows)


with OUT_MEMBERS.open("w", newline="") as f:

    fields = [
        "orthogroup",
        "cluster_index",
        "sequence_id",
        "accession",
        "protein_id",
        "original_tmd_index",
        "alignment_start",
        "alignment_end",
        "native_start",
        "native_end",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(member_rows)


# --------------------------------------------------
# Summary
# --------------------------------------------------

high_support = [
    r for r in cluster_rows
    if r["support_fraction"] >= 0.8
]

clean_high_support = [
    r for r in high_support
    if r["proteins_with_multiple_tmds_in_cluster"] == 0
]

by_family = defaultdict(int)

for r in clean_high_support:
    by_family[r["orthogroup"]] += 1

print("Families processed:", len(families))
print("Total overlap clusters:", len(cluster_rows))
print("Clusters with support >=80%:", len(high_support))
print(
    "Clean clusters with support >=80%:",
    len(clean_high_support)
)
print(
    "Families with >=3 clean high-support clusters:",
    sum(n >= 3 for n in by_family.values())
)

print("Output clusters:", OUT_CLUSTERS)
print("Output members:", OUT_MEMBERS)
