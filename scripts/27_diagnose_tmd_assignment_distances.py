from pathlib import Path
from collections import defaultdict
import csv
import statistics

ASSIGN = Path(
    "work/codon/tmd_assignment/segment_tmd_assignments.tsv"
)

UNITS = Path(
    "work/codon/tmd_assignment/homologous_tmd_unit_summary.tsv"
)

rows = []

with ASSIGN.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        r["absolute_anchor_distance"] = float(
            r["absolute_anchor_distance"]
        )
        r["anchor_relative_position"] = float(
            r["anchor_relative_position"]
        )
        rows.append(r)


print("Assignment-distance diagnostic")
print()

for anchor in ["start", "end", "center"]:

    x = [
        r["absolute_anchor_distance"]
        for r in rows
        if r["anchor_type"] == anchor
    ]

    xs = sorted(x)

    def q(p):
        if not xs:
            return None
        i = int(round((len(xs) - 1) * p))
        return xs[i]

    print("Anchor:", anchor)
    print("  assignments:", len(xs))
    print("  median abs distance:", statistics.median(xs))
    print("  75%:", q(0.75))
    print("  90%:", q(0.90))
    print("  95%:", q(0.95))
    print("  99%:", q(0.99))
    print("  max:", max(xs))

    for cutoff in [25, 50, 100, 200]:
        n = sum(v > cutoff for v in xs)
        print(
            f"  >{cutoff} aa:",
            n,
            f"({n/len(xs):.3%})"
        )

    print()


# -----------------------------------------------
# qualifying units with largest variance
# -----------------------------------------------

units = []

with UNITS.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for r in reader:

        if r["qualifying_unit"] != "1":
            continue

        if r["sample_variance"] == "":
            continue

        r["sample_variance"] = float(
            r["sample_variance"]
        )

        units.append(r)


print("Top qualifying units by variance")
print()

for anchor in ["start", "end", "center"]:

    sub = [
        r for r in units
        if r["anchor_type"] == anchor
    ]

    sub.sort(
        key=lambda r: r["sample_variance"],
        reverse=True
    )

    print("Anchor:", anchor)

    for r in sub[:10]:
        print(
            r["orthogroup"],
            "cluster", r["cluster_index"],
            "positive_species=",
            r["segment_positive_species_n"],
            "variance=",
            f'{r["sample_variance"]:.2f}'
        )

    print()


# -----------------------------------------------
# show member assignments for top 5 start units
# -----------------------------------------------

top_start = sorted(
    [
        r for r in units
        if r["anchor_type"] == "start"
    ],
    key=lambda r: r["sample_variance"],
    reverse=True
)[:5]

print("Detailed assignments for top 5 START units")
print()

for unit in top_start:

    og = unit["orthogroup"]
    ci = unit["cluster_index"]

    print(
        og,
        "cluster", ci,
        "variance=",
        unit["sample_variance"]
    )

    members = [
        r for r in rows
        if r["anchor_type"] == "start"
        and r["orthogroup"] == og
        and r["cluster_index"] == ci
    ]

    members.sort(
        key=lambda r: r["accession"]
    )

    for r in members:

        print(
            " ",
            r["accession"],
            r["protein_id"],
            "segment_center=",
            r["segment_center"],
            "TMD_start=",
            r["native_tmd_start"],
            "relative=",
            r["anchor_relative_position"],
            "abs=",
            r["absolute_anchor_distance"],
        )

    print()
