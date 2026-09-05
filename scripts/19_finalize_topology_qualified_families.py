from pathlib import Path
from collections import defaultdict
import csv

CLUSTERS = Path(
    "work/topology/merged/tmd_overlap_clusters.tsv"
)

MEMBERS = Path(
    "work/topology/merged/tmd_overlap_cluster_members.tsv"
)

OUT_FAMILIES = Path(
    "work/topology/merged/topology_qualified_families.tsv"
)

OUT_TMDS = Path(
    "work/topology/merged/topology_qualified_tmd_clusters.tsv"
)

OUT_EXCLUDED = Path(
    "work/topology/merged/topology_excluded_families.tsv"
)

MIN_SUPPORT = 0.80
MAX_START_MAD = 2.0
MAX_END_MAD = 2.0
MIN_GOOD_CLUSTERS = 3


# --------------------------------------------------
# 1. Load clusters
# --------------------------------------------------

clusters = []
clusters_by_family = defaultdict(list)

with CLUSTERS.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        row["cluster_index"] = int(row["cluster_index"])
        row["family_member_count"] = int(row["family_member_count"])
        row["support_n"] = int(row["support_n"])
        row["support_fraction"] = float(row["support_fraction"])
        row["proteins_with_multiple_tmds_in_cluster"] = int(
            row["proteins_with_multiple_tmds_in_cluster"]
        )
        row["start_mad"] = float(row["start_mad"])
        row["end_mad"] = float(row["end_mad"])
        row["median_alignment_start"] = float(
            row["median_alignment_start"]
        )
        row["median_alignment_end"] = float(
            row["median_alignment_end"]
        )

        row["passes_cluster_qc"] = (
            row["support_fraction"] >= MIN_SUPPORT
            and row["proteins_with_multiple_tmds_in_cluster"] == 0
            and row["start_mad"] <= MAX_START_MAD
            and row["end_mad"] <= MAX_END_MAD
        )

        clusters.append(row)
        clusters_by_family[row["orthogroup"]].append(row)


# --------------------------------------------------
# 2. Load member assignments for order checking
# --------------------------------------------------

members_by_family_seq = defaultdict(lambda: defaultdict(list))

with MEMBERS.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        og = row["orthogroup"]
        seq = row["sequence_id"]

        members_by_family_seq[og][seq].append({
            "cluster_index": int(row["cluster_index"]),
            "original_tmd_index": int(row["original_tmd_index"]),
        })


# --------------------------------------------------
# 3. Evaluate family order consistency
# --------------------------------------------------

family_rows = []
excluded_rows = []
qualified_clusters = []

for og, fam_clusters in clusters_by_family.items():

    good = [
        x for x in fam_clusters
        if x["passes_cluster_qc"]
    ]

    good_cluster_ids = {
        x["cluster_index"] for x in good
    }

    order_violations = 0
    sequences_checked = 0

    for seq_id, assignments in members_by_family_seq[og].items():

        x = [
            a for a in assignments
            if a["cluster_index"] in good_cluster_ids
        ]

        if len(x) < 2:
            continue

        sequences_checked += 1

        # Sort according to native N -> C TMD order
        x.sort(key=lambda a: a["original_tmd_index"])

        cluster_order = [
            a["cluster_index"] for a in x
        ]

        # Because clusters were numbered by median
        # alignment position, this should increase.
        if any(
            cluster_order[i] >= cluster_order[i + 1]
            for i in range(len(cluster_order) - 1)
        ):
            order_violations += 1

    qualifies = (
        len(good) >= MIN_GOOD_CLUSTERS
        and order_violations == 0
    )

    row = {
        "orthogroup": og,
        "total_clusters": len(fam_clusters),
        "qc_pass_clusters": len(good),
        "sequences_checked_for_order": sequences_checked,
        "order_violations": order_violations,
        "topology_qualified": int(qualifies),
    }

    family_rows.append(row)

    if qualifies:
        for c in good:
            qualified_clusters.append(c)
    else:
        reasons = []

        if len(good) < MIN_GOOD_CLUSTERS:
            reasons.append("fewer_than_3_qc_pass_clusters")

        if order_violations > 0:
            reasons.append("tmd_order_violation")

        excluded_rows.append({
            **row,
            "exclusion_reason": ";".join(reasons),
        })


# --------------------------------------------------
# 4. Sort
# --------------------------------------------------

def og_num(row):
    return int(row["orthogroup"].split("_")[-1])


family_rows.sort(key=og_num)
excluded_rows.sort(key=og_num)

qualified_clusters.sort(
    key=lambda r: (
        int(r["orthogroup"].split("_")[-1]),
        r["cluster_index"],
    )
)


# --------------------------------------------------
# 5. Write family table
# --------------------------------------------------

family_fields = [
    "orthogroup",
    "total_clusters",
    "qc_pass_clusters",
    "sequences_checked_for_order",
    "order_violations",
    "topology_qualified",
]

with OUT_FAMILIES.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=family_fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(family_rows)


# --------------------------------------------------
# 6. Write qualified TMD clusters
# --------------------------------------------------

cluster_fields = [
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

with OUT_TMDS.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=cluster_fields,
        delimiter="\t",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(qualified_clusters)


# --------------------------------------------------
# 7. Write excluded families
# --------------------------------------------------

excluded_fields = family_fields + ["exclusion_reason"]

with OUT_EXCLUDED.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=excluded_fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(excluded_rows)


# --------------------------------------------------
# 8. Summary
# --------------------------------------------------

qualified = [
    r for r in family_rows
    if r["topology_qualified"] == 1
]

print("Families evaluated:", len(family_rows))
print("Topology-qualified families:", len(qualified))
print("Excluded families:", len(excluded_rows))
print(
    "Qualified homologous TMD clusters:",
    len(qualified_clusters)
)
print(
    "Families with any TMD-order violation:",
    sum(r["order_violations"] > 0 for r in family_rows)
)

print()
print("Outputs:")
print(OUT_FAMILIES)
print(OUT_TMDS)
print(OUT_EXCLUDED)

if excluded_rows:
    print()
    print("Excluded families:")
    for r in excluded_rows:
        print(
            r["orthogroup"],
            "good_clusters=" + str(r["qc_pass_clusters"]),
            "order_violations=" + str(r["order_violations"]),
            r["exclusion_reason"],
        )
