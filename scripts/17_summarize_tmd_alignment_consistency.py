from pathlib import Path
from collections import defaultdict
import csv
import statistics

INFILE = Path(
    "work/topology/merged/tmd_alignment_coordinates.tsv"
)

FAMILY_SUMMARY = Path(
    "work/topology/merged/orthogroup_topology_summary.tsv"
)

OUTFILE = Path(
    "work/topology/merged/tmd_alignment_consistency.tsv"
)

# --------------------------------------------------
# Load family member counts
# --------------------------------------------------

member_counts = {}

with FAMILY_SUMMARY.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        member_counts[row["orthogroup"]] = int(row["member_count"])

# --------------------------------------------------
# Load TMD alignment coordinates
# --------------------------------------------------

groups = defaultdict(list)

with INFILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        row["tmd_index"] = int(row["tmd_index"])
        row["alignment_start"] = int(row["alignment_start"])
        row["alignment_end"] = int(row["alignment_end"])
        row["native_length"] = int(row["native_length"])

        key = (row["orthogroup"], row["tmd_index"])
        groups[key].append(row)


def median_abs_deviation(values):
    med = statistics.median(values)
    return statistics.median(
        [abs(x - med) for x in values]
    )


rows = []

for (og, tmd_index), records in groups.items():

    starts = [x["alignment_start"] for x in records]
    ends = [x["alignment_end"] for x in records]
    lengths = [x["native_length"] for x in records]

    support_n = len(records)
    family_n = member_counts[og]

    rows.append({
        "orthogroup": og,
        "tmd_index": tmd_index,
        "family_member_count": family_n,
        "support_n": support_n,
        "support_fraction": support_n / family_n,

        "median_alignment_start": statistics.median(starts),
        "min_alignment_start": min(starts),
        "max_alignment_start": max(starts),
        "start_range": max(starts) - min(starts),
        "start_mad": median_abs_deviation(starts),

        "median_alignment_end": statistics.median(ends),
        "min_alignment_end": min(ends),
        "max_alignment_end": max(ends),
        "end_range": max(ends) - min(ends),
        "end_mad": median_abs_deviation(ends),

        "median_tmd_length": statistics.median(lengths),
        "min_tmd_length": min(lengths),
        "max_tmd_length": max(lengths),
    })

rows.sort(
    key=lambda x: (
        int(x["orthogroup"].split("_")[-1]),
        x["tmd_index"]
    )
)

with OUTFILE.open("w", newline="") as f:

    fieldnames = [
        "orthogroup",
        "tmd_index",
        "family_member_count",
        "support_n",
        "support_fraction",

        "median_alignment_start",
        "min_alignment_start",
        "max_alignment_start",
        "start_range",
        "start_mad",

        "median_alignment_end",
        "min_alignment_end",
        "max_alignment_end",
        "end_range",
        "end_mad",

        "median_tmd_length",
        "min_tmd_length",
        "max_tmd_length",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

print("TMD-index units summarized:", len(rows))
print("Families represented:", len(set(r["orthogroup"] for r in rows)))

print(
    "Units with full family support:",
    sum(r["support_n"] == r["family_member_count"] for r in rows)
)

print(
    "Units with support >=80%:",
    sum(r["support_fraction"] >= 0.8 for r in rows)
)

print("Output:", OUTFILE)
