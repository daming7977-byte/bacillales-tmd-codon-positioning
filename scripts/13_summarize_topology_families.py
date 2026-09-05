from pathlib import Path
from collections import defaultdict, Counter
import csv
import statistics

INFILE = Path(
    "work/topology/merged/deeptmhmm_all_predictions.tsv"
)

OUTFILE = Path(
    "work/topology/merged/orthogroup_topology_summary.tsv"
)

groups = defaultdict(list)

with INFILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        row["n_tmd"] = int(row["n_tmd"])
        groups[row["orthogroup"]].append(row)

rows = []

for og, members in groups.items():

    n_tmds = [x["n_tmd"] for x in members]

    counts = Counter(n_tmds)
    modal_n_tmd = sorted(
        counts.items(),
        key=lambda x: (-x[1], x[0])
    )[0][0]

    n_with_tmd = sum(x >= 1 for x in n_tmds)
    n_ge3_tmd = sum(x >= 3 for x in n_tmds)

    rows.append({
        "orthogroup": og,
        "member_count": len(members),
        "n_with_tmd": n_with_tmd,
        "n_ge3_tmd": n_ge3_tmd,
        "fraction_with_tmd": n_with_tmd / len(members),
        "fraction_ge3_tmd": n_ge3_tmd / len(members),
        "min_n_tmd": min(n_tmds),
        "median_n_tmd": statistics.median(n_tmds),
        "max_n_tmd": max(n_tmds),
        "modal_n_tmd": modal_n_tmd,
        "n_tmd_distribution": ";".join(
            f"{k}:{v}" for k, v in sorted(counts.items())
        ),
    })

def og_number(x):
    try:
        return int(x["orthogroup"].split("_")[-1])
    except Exception:
        return 10**12

rows.sort(key=og_number)

with OUTFILE.open("w", newline="") as f:
    fieldnames = [
        "orthogroup",
        "member_count",
        "n_with_tmd",
        "n_ge3_tmd",
        "fraction_with_tmd",
        "fraction_ge3_tmd",
        "min_n_tmd",
        "median_n_tmd",
        "max_n_tmd",
        "modal_n_tmd",
        "n_tmd_distribution",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

print("Orthogroups:", len(rows))

print(
    "Families with >=1 member having >=3 TMD:",
    sum(r["n_ge3_tmd"] >= 1 for r in rows)
)

print(
    "Families with >=50% members having >=3 TMD:",
    sum(r["fraction_ge3_tmd"] >= 0.5 for r in rows)
)

print(
    "Families with all members having >=3 TMD:",
    sum(r["n_ge3_tmd"] == r["member_count"] for r in rows)
)

print(
    "Families with modal nTMD >=3:",
    sum(r["modal_n_tmd"] >= 3 for r in rows)
)

print("Output:", OUTFILE)
