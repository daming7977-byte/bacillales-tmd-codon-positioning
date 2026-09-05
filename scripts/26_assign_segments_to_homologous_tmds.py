from pathlib import Path
from collections import defaultdict, Counter
import csv
import statistics

SEGMENT_FILE = Path(
    "work/codon/segments/low_adaptation_segments.tsv"
)

QUALIFIED_CLUSTER_FILE = Path(
    "work/topology/merged/topology_qualified_tmd_clusters.tsv"
)

CLUSTER_MEMBER_FILE = Path(
    "work/topology/merged/tmd_overlap_cluster_members.tsv"
)

OUTDIR = Path("work/codon/tmd_assignment")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_ASSIGN = OUTDIR / "segment_tmd_assignments.tsv"
OUT_UNITS = OUTDIR / "homologous_tmd_unit_summary.tsv"


# --------------------------------------------------
# 1. Load qualified homologous TMD clusters
# --------------------------------------------------

qualified_clusters = {}
qualified_keys = set()

with QUALIFIED_CLUSTER_FILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        og = row["orthogroup"]
        cluster = int(row["cluster_index"])
        key = (og, cluster)

        qualified_keys.add(key)

        qualified_clusters[key] = {
            "orthogroup": og,
            "cluster_index": cluster,
            "family_member_count": int(row["family_member_count"]),
            "support_n": int(row["support_n"]),
            "support_fraction": float(row["support_fraction"]),
            "start_mad": float(row["start_mad"]),
            "end_mad": float(row["end_mad"]),
        }

print(
    "Qualified homologous TMD clusters:",
    len(qualified_clusters)
)


# --------------------------------------------------
# 2. Load member-specific native TMD coordinates
#    but only for the 1086 frozen qualified clusters
# --------------------------------------------------

members_by_sequence = defaultdict(list)
member_count = 0

with CLUSTER_MEMBER_FILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        og = row["orthogroup"]
        cluster = int(row["cluster_index"])
        key = (og, cluster)

        if key not in qualified_keys:
            continue

        start = int(row["native_start"])
        end = int(row["native_end"])

        rec = {
            "orthogroup": og,
            "cluster_index": cluster,
            "sequence_id": row["sequence_id"],
            "accession": row["accession"],
            "protein_id": row["protein_id"],
            "native_start": start,
            "native_end": end,
            "native_center": (start + end) / 2,
            "original_tmd_index": int(row["original_tmd_index"]),
        }

        members_by_sequence[row["sequence_id"]].append(rec)
        member_count += 1

print(
    "Qualified cluster member records:",
    member_count
)
print(
    "Proteins with >=1 qualified TMD:",
    len(members_by_sequence)
)


# --------------------------------------------------
# 3. Load observed low-adaptation segments
# --------------------------------------------------

segments_by_sequence = defaultdict(list)
segment_total = 0

with SEGMENT_FILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        rec = {
            "orthogroup": row["orthogroup"],
            "accession": row["accession"],
            "protein_id": row["protein_id"],
            "sequence_id": row["sequence_id"],
            "segment_index": int(row["segment_index"]),
            "segment_start": int(row["segment_start"]),
            "segment_end": int(row["segment_end"]),
            "segment_center": float(row["segment_center"]),
            "segment_length_codons": int(
                row["segment_length_codons"]
            ),
        }

        segments_by_sequence[row["sequence_id"]].append(rec)
        segment_total += 1

print("Observed low-adaptation segments:", segment_total)


# --------------------------------------------------
# 4. Independent assignment for start/end/center
#
#    For every segment:
#      find nearest qualified TMD anchor
#
#    Then:
#      if multiple segments from same protein are
#      assigned to same cluster, retain nearest one.
# --------------------------------------------------

anchor_fields = {
    "start": "native_start",
    "end": "native_end",
    "center": "native_center",
}

assignment_rows = []

# retained[(anchor_type, sequence_id, og, cluster)]
#          = best assignment
retained = {}

segments_without_anchor = 0

for sequence_id, segments in segments_by_sequence.items():

    tmds = members_by_sequence.get(sequence_id, [])

    if not tmds:
        segments_without_anchor += len(segments)
        continue

    for seg in segments:

        for anchor_type, anchor_field in anchor_fields.items():

            # Find closest anchor.
            # Deterministic tie-break:
            # lower cluster_index wins.
            candidates = []

            for tmd in tmds:

                anchor_pos = tmd[anchor_field]

                relative = (
                    seg["segment_center"] - anchor_pos
                )

                candidates.append(
                    (
                        abs(relative),
                        tmd["cluster_index"],
                        relative,
                        anchor_pos,
                        tmd,
                    )
                )

            candidates.sort(
                key=lambda x: (x[0], x[1])
            )

            (
                abs_distance,
                cluster_index,
                relative,
                anchor_pos,
                tmd,
            ) = candidates[0]

            key = (
                anchor_type,
                sequence_id,
                tmd["orthogroup"],
                cluster_index,
            )

            candidate = {
                "anchor_type": anchor_type,
                "orthogroup": tmd["orthogroup"],
                "cluster_index": cluster_index,
                "sequence_id": sequence_id,
                "accession": seg["accession"],
                "protein_id": seg["protein_id"],
                "segment_index": seg["segment_index"],
                "segment_start": seg["segment_start"],
                "segment_end": seg["segment_end"],
                "segment_center": seg["segment_center"],
                "segment_length_codons":
                    seg["segment_length_codons"],
                "native_tmd_start": tmd["native_start"],
                "native_tmd_end": tmd["native_end"],
                "native_tmd_center": tmd["native_center"],
                "anchor_position": anchor_pos,
                "anchor_relative_position": relative,
                "absolute_anchor_distance": abs_distance,
                "original_tmd_index":
                    tmd["original_tmd_index"],
            }

            # Same protein + same TMD cluster:
            # retain closest segment only.
            if key not in retained:
                retained[key] = candidate
            else:
                old = retained[key]

                new_key = (
                    candidate["absolute_anchor_distance"],
                    candidate["segment_index"],
                )

                old_key = (
                    old["absolute_anchor_distance"],
                    old["segment_index"],
                )

                if new_key < old_key:
                    retained[key] = candidate


assignment_rows = list(retained.values())

assignment_rows.sort(
    key=lambda r: (
        r["anchor_type"],
        r["orthogroup"],
        r["cluster_index"],
        r["accession"],
        r["sequence_id"],
    )
)

print(
    "Segments lacking any qualified TMD anchor:",
    segments_without_anchor
)

print(
    "Retained protein-TMD assignments:",
    len(assignment_rows)
)


# --------------------------------------------------
# 5. Write assignment table
# --------------------------------------------------

with OUT_ASSIGN.open("w", newline="") as f:

    fields = [
        "anchor_type",
        "orthogroup",
        "cluster_index",
        "sequence_id",
        "accession",
        "protein_id",
        "segment_index",
        "segment_start",
        "segment_end",
        "segment_center",
        "segment_length_codons",
        "native_tmd_start",
        "native_tmd_end",
        "native_tmd_center",
        "anchor_position",
        "anchor_relative_position",
        "absolute_anchor_distance",
        "original_tmd_index",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(assignment_rows)


# --------------------------------------------------
# 6. Summarize homologous TMD units
# --------------------------------------------------

assignments_by_unit = defaultdict(list)

for row in assignment_rows:
    key = (
        row["anchor_type"],
        row["orthogroup"],
        row["cluster_index"],
    )

    assignments_by_unit[key].append(row)


unit_rows = []

for anchor_type in anchor_fields:

    for key in sorted(qualified_clusters):

        og, cluster = key
        qc = qualified_clusters[key]

        assigned = assignments_by_unit.get(
            (anchor_type, og, cluster),
            []
        )

        positive_species = sorted({
            x["accession"]
            for x in assigned
        })

        relative_positions = [
            float(x["anchor_relative_position"])
            for x in assigned
        ]

        support_n = qc["support_n"]

        represented_ge8 = int(
            support_n >= 8
        )

        segment_species_n = len(positive_species)

        qualifying = int(
            support_n >= 8
            and segment_species_n >= 3
        )

        if segment_species_n >= 2:
            sample_variance = statistics.variance(
                relative_positions
            )
        else:
            sample_variance = ""

        if relative_positions:
            median_relative_position = statistics.median(
                relative_positions
            )
        else:
            median_relative_position = ""

        unit_rows.append({
            "anchor_type": anchor_type,
            "orthogroup": og,
            "cluster_index": cluster,
            "family_member_count":
                qc["family_member_count"],
            "tmd_species_support_n": support_n,
            "tmd_support_fraction":
                qc["support_fraction"],
            "start_mad": qc["start_mad"],
            "end_mad": qc["end_mad"],
            "represented_ge8_species":
                represented_ge8,
            "segment_positive_species_n":
                segment_species_n,
            "segment_positive_species":
                ";".join(positive_species),
            "median_anchor_relative_position":
                median_relative_position,
            "sample_variance":
                sample_variance,
            "qualifying_unit": qualifying,
        })


with OUT_UNITS.open("w", newline="") as f:

    fields = [
        "anchor_type",
        "orthogroup",
        "cluster_index",
        "family_member_count",
        "tmd_species_support_n",
        "tmd_support_fraction",
        "start_mad",
        "end_mad",
        "represented_ge8_species",
        "segment_positive_species_n",
        "segment_positive_species",
        "median_anchor_relative_position",
        "sample_variance",
        "qualifying_unit",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(unit_rows)


# --------------------------------------------------
# 7. Main observed summaries
# --------------------------------------------------

print()
print("Observed homologous-TMD summary")

for anchor_type in ["start", "end", "center"]:

    rows = [
        r for r in unit_rows
        if r["anchor_type"] == anchor_type
    ]

    represented = [
        r for r in rows
        if r["represented_ge8_species"] == 1
    ]

    qualifying = [
        r for r in rows
        if r["qualifying_unit"] == 1
    ]

    variances = [
        float(r["sample_variance"])
        for r in qualifying
        if r["sample_variance"] != ""
    ]

    families = {
        r["orthogroup"]
        for r in qualifying
    }

    print()
    print("Anchor:", anchor_type)
    print(
        "  TMD units represented in >=8 species:",
        len(represented)
    )
    print(
        "  Qualifying units (>=3 segment-positive species):",
        len(qualifying)
    )
    print(
        "  Orthogroups contributing qualifying units:",
        len(families)
    )

    if variances:
        print(
            "  Median sample variance:",
            f"{statistics.median(variances):.6f}"
        )

        print(
            "  Min sample variance:",
            f"{min(variances):.6f}"
        )

        print(
            "  Max sample variance:",
            f"{max(variances):.6f}"
        )


# --------------------------------------------------
# 8. Distribution of segment-positive species
# --------------------------------------------------

print()
print(
    "Distribution of segment-positive species "
    "per TMD unit"
)

for anchor_type in ["start", "end", "center"]:

    counts = Counter()

    for r in unit_rows:

        if r["anchor_type"] != anchor_type:
            continue

        if not r["represented_ge8_species"]:
            continue

        counts[r["segment_positive_species_n"]] += 1

    print()
    print(anchor_type)

    for n in sorted(counts):
        print(
            f"  {n} species: {counts[n]}"
        )


print()
print("Outputs:")
print(OUT_ASSIGN)
print(OUT_UNITS)
