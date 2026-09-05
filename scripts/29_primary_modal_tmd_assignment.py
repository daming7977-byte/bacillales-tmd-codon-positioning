from pathlib import Path
from collections import defaultdict, Counter
import csv
import statistics

SEGMENT_FILE = Path(
    "work/codon/segments/low_adaptation_segments.tsv"
)

MODAL_FILE = Path(
    "work/topology/merged/modal_tmd_primary_set_diagnostic.tsv"
)

ALL_CLUSTER_MEMBERS = Path(
    "work/topology/merged/tmd_overlap_cluster_members.tsv"
)

QUALIFIED_CLUSTERS = Path(
    "work/topology/merged/topology_qualified_tmd_clusters.tsv"
)

OUTDIR = Path("work/codon/primary_tmd_assignment")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_ASSIGN_ALL = OUTDIR / "all_nearest_tmd_assignments.tsv"
OUT_ASSIGN_PRIMARY = OUTDIR / "primary_qualified_tmd_assignments.tsv"
OUT_UNITS = OUTDIR / "primary_homologous_tmd_unit_summary.tsv"


# --------------------------------------------------
# 1. Modal-matching primary proteins
# --------------------------------------------------

modal_ok = set()
modal_all = set()

with MODAL_FILE.open() as f:
    r = csv.DictReader(f, delimiter="\t")

    for x in r:
        sid = x["sequence_id"]
        modal_all.add(sid)

        if x["matches_modal"] == "1":
            modal_ok.add(sid)

print("Qualified-family proteins:", len(modal_all))
print("Modal-matching proteins:", len(modal_ok))


# --------------------------------------------------
# 2. Frozen qualified homologous clusters
# --------------------------------------------------

qualified_keys = set()
qualified_qc = {}

with QUALIFIED_CLUSTERS.open() as f:
    r = csv.DictReader(f, delimiter="\t")

    for x in r:
        key = (
            x["orthogroup"],
            int(x["cluster_index"])
        )

        qualified_keys.add(key)

        qualified_qc[key] = {
            "start_mad": float(x["start_mad"]),
            "end_mad": float(x["end_mad"]),
        }

print(
    "Frozen qualified homologous TMD clusters:",
    len(qualified_keys)
)


# --------------------------------------------------
# 3. ALL TMDs in modal-matching proteins
#
# Important:
# nearest-anchor assignment is performed BEFORE
# qualified-cluster filtering.
# --------------------------------------------------

tmds_by_sequence = defaultdict(list)

# filtered support of each homologous cluster
modal_cluster_species = defaultdict(set)

all_member_records = 0
modal_member_records = 0

with ALL_CLUSTER_MEMBERS.open() as f:
    r = csv.DictReader(f, delimiter="\t")

    for x in r:

        all_member_records += 1

        sid = x["sequence_id"]

        if sid not in modal_ok:
            continue

        modal_member_records += 1

        start = int(x["native_start"])
        end = int(x["native_end"])

        rec = {
            "orthogroup": x["orthogroup"],
            "cluster_index": int(x["cluster_index"]),
            "sequence_id": sid,
            "accession": x["accession"],
            "protein_id": x["protein_id"],
            "original_tmd_index":
                int(x["original_tmd_index"]),
            "native_start": start,
            "native_end": end,
            "native_center": (start + end) / 2,
        }

        tmds_by_sequence[sid].append(rec)

        key = (
            x["orthogroup"],
            int(x["cluster_index"])
        )

        modal_cluster_species[key].add(
            x["accession"]
        )


print("All cluster-member records:", all_member_records)
print(
    "Modal-matching cluster-member records:",
    modal_member_records
)
print(
    "Modal proteins with >=1 TMD:",
    len(tmds_by_sequence)
)


# --------------------------------------------------
# 4. Load segments but keep only modal proteins
# --------------------------------------------------

segments_by_sequence = defaultdict(list)

segments_total = 0
segments_modal = 0

with SEGMENT_FILE.open() as f:
    r = csv.DictReader(f, delimiter="\t")

    for x in r:

        segments_total += 1

        sid = x["sequence_id"]

        if sid not in modal_ok:
            continue

        segments_modal += 1

        segments_by_sequence[sid].append({
            "orthogroup": x["orthogroup"],
            "accession": x["accession"],
            "protein_id": x["protein_id"],
            "sequence_id": sid,
            "segment_index":
                int(x["segment_index"]),
            "segment_start":
                int(x["segment_start"]),
            "segment_end":
                int(x["segment_end"]),
            "segment_center":
                float(x["segment_center"]),
            "segment_length_codons":
                int(x["segment_length_codons"]),
        })

print("All observed segments:", segments_total)
print(
    "Segments in modal-matching proteins:",
    segments_modal
)


# --------------------------------------------------
# 5. Assign each segment to nearest ALL-TMD anchor
# --------------------------------------------------

anchor_fields = {
    "start": "native_start",
    "end": "native_end",
    "center": "native_center",
}

# Best segment retained per
# anchor + protein + assigned TMD cluster
retained = {}

no_tmd = 0

for sid, segments in segments_by_sequence.items():

    tmds = tmds_by_sequence.get(sid, [])

    if not tmds:
        no_tmd += len(segments)
        continue

    for seg in segments:

        for anchor_type, field in anchor_fields.items():

            candidates = []

            for tmd in tmds:

                anchor = tmd[field]

                rel = (
                    seg["segment_center"] - anchor
                )

                candidates.append(
                    (
                        abs(rel),
                        tmd["original_tmd_index"],
                        tmd["cluster_index"],
                        rel,
                        anchor,
                        tmd,
                    )
                )

            # deterministic:
            # nearest distance,
            # then lower native TMD index,
            # then lower cluster index
            candidates.sort(
                key=lambda z: (
                    z[0],
                    z[1],
                    z[2],
                )
            )

            (
                absdist,
                original_index,
                cluster_index,
                rel,
                anchor,
                tmd,
            ) = candidates[0]

            key = (
                anchor_type,
                sid,
                tmd["orthogroup"],
                cluster_index,
            )

            row = {
                "anchor_type": anchor_type,
                "orthogroup": tmd["orthogroup"],
                "cluster_index": cluster_index,
                "cluster_is_qualified":
                    int(
                        (
                            tmd["orthogroup"],
                            cluster_index
                        )
                        in qualified_keys
                    ),
                "sequence_id": sid,
                "accession": seg["accession"],
                "protein_id": seg["protein_id"],
                "segment_index":
                    seg["segment_index"],
                "segment_start":
                    seg["segment_start"],
                "segment_end":
                    seg["segment_end"],
                "segment_center":
                    seg["segment_center"],
                "segment_length_codons":
                    seg["segment_length_codons"],
                "original_tmd_index":
                    original_index,
                "native_tmd_start":
                    tmd["native_start"],
                "native_tmd_end":
                    tmd["native_end"],
                "native_tmd_center":
                    tmd["native_center"],
                "anchor_position": anchor,
                "anchor_relative_position": rel,
                "absolute_anchor_distance":
                    absdist,
            }

            if key not in retained:
                retained[key] = row
            else:
                old = retained[key]

                if (
                    row["absolute_anchor_distance"],
                    row["segment_index"],
                ) < (
                    old["absolute_anchor_distance"],
                    old["segment_index"],
                ):
                    retained[key] = row


all_assignments = list(retained.values())

all_assignments.sort(
    key=lambda x: (
        x["anchor_type"],
        x["orthogroup"],
        x["cluster_index"],
        x["accession"],
    )
)

print(
    "Segments without any TMD:",
    no_tmd
)
print(
    "Retained nearest-TMD assignments:",
    len(all_assignments)
)


# --------------------------------------------------
# 6. Write all nearest assignments
# --------------------------------------------------

fields_assign = [
    "anchor_type",
    "orthogroup",
    "cluster_index",
    "cluster_is_qualified",
    "sequence_id",
    "accession",
    "protein_id",
    "segment_index",
    "segment_start",
    "segment_end",
    "segment_center",
    "segment_length_codons",
    "original_tmd_index",
    "native_tmd_start",
    "native_tmd_end",
    "native_tmd_center",
    "anchor_position",
    "anchor_relative_position",
    "absolute_anchor_distance",
]

with OUT_ASSIGN_ALL.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields_assign,
        delimiter="\t"
    )
    w.writeheader()
    w.writerows(all_assignments)


# --------------------------------------------------
# 7. Primary assignments:
#    nearest TMD itself must belong to frozen
#    qualified homologous cluster.
#
# No reassignment if nearest cluster failed QC.
# --------------------------------------------------

primary_assignments = [
    x for x in all_assignments
    if x["cluster_is_qualified"] == 1
]

with OUT_ASSIGN_PRIMARY.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields_assign,
        delimiter="\t"
    )
    w.writeheader()
    w.writerows(primary_assignments)


print(
    "Assignments whose nearest TMD is qualified:",
    len(primary_assignments)
)

print(
    "Assignments excluded because nearest TMD "
    "cluster is non-qualified:",
    len(all_assignments)
    - len(primary_assignments)
)


# --------------------------------------------------
# 8. Unit-level primary analysis
#
# Recompute species support AFTER modal filter.
# --------------------------------------------------

assign_by_unit = defaultdict(list)

for x in primary_assignments:

    key = (
        x["anchor_type"],
        x["orthogroup"],
        x["cluster_index"],
    )

    assign_by_unit[key].append(x)


unit_rows = []

for anchor_type in ["start", "end", "center"]:

    for key in sorted(qualified_keys):

        og, ci = key

        modal_species = sorted(
            modal_cluster_species.get(
                key, set()
            )
        )

        support_n = len(modal_species)

        assigned = assign_by_unit.get(
            (anchor_type, og, ci),
            []
        )

        positive_species = sorted({
            x["accession"]
            for x in assigned
        })

        pos = [
            float(
                x["anchor_relative_position"]
            )
            for x in assigned
        ]

        qualifying = int(
            support_n >= 8
            and len(positive_species) >= 3
        )

        variance = (
            statistics.variance(pos)
            if len(pos) >= 2
            else ""
        )

        medpos = (
            statistics.median(pos)
            if pos
            else ""
        )

        unit_rows.append({
            "anchor_type": anchor_type,
            "orthogroup": og,
            "cluster_index": ci,
            "modal_species_support_n":
                support_n,
            "modal_species":
                ";".join(modal_species),
            "represented_ge8_species":
                int(support_n >= 8),
            "segment_positive_species_n":
                len(positive_species),
            "segment_positive_species":
                ";".join(positive_species),
            "median_anchor_relative_position":
                medpos,
            "sample_variance":
                variance,
            "qualifying_unit":
                qualifying,
        })


fields_units = [
    "anchor_type",
    "orthogroup",
    "cluster_index",
    "modal_species_support_n",
    "modal_species",
    "represented_ge8_species",
    "segment_positive_species_n",
    "segment_positive_species",
    "median_anchor_relative_position",
    "sample_variance",
    "qualifying_unit",
]

with OUT_UNITS.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields_units,
        delimiter="\t"
    )
    w.writeheader()
    w.writerows(unit_rows)


# --------------------------------------------------
# 9. Diagnostics
# --------------------------------------------------

print()
print("CORRECTED PRIMARY OBSERVED ANALYSIS")

for anchor in ["start", "end", "center"]:

    u = [
        x for x in unit_rows
        if x["anchor_type"] == anchor
    ]

    represented = [
        x for x in u
        if x["represented_ge8_species"] == 1
    ]

    qualifying = [
        x for x in u
        if x["qualifying_unit"] == 1
    ]

    variances = [
        float(x["sample_variance"])
        for x in qualifying
        if x["sample_variance"] != ""
    ]

    fams = {
        x["orthogroup"]
        for x in qualifying
    }

    arows = [
        x for x in primary_assignments
        if x["anchor_type"] == anchor
    ]

    distances = sorted(
        float(x["absolute_anchor_distance"])
        for x in arows
    )

    print()
    print("Anchor:", anchor)
    print(
        "  qualified units with >=8 modal species:",
        len(represented)
    )
    print(
        "  qualifying units (>=3 positive species):",
        len(qualifying)
    )
    print(
        "  contributing orthogroups:",
        len(fams)
    )

    if variances:
        print(
            "  median sample variance:",
            f"{statistics.median(variances):.6f}"
        )

    if distances:
        print(
            "  assignment median abs distance:",
            statistics.median(distances)
        )

        print(
            "  assignments >50 aa:",
            sum(x > 50 for x in distances)
        )

        print(
            "  assignments >100 aa:",
            sum(x > 100 for x in distances)
        )


print()
print("Top 10 START qualifying units by variance")

top = [
    x for x in unit_rows
    if x["anchor_type"] == "start"
    and x["qualifying_unit"] == 1
]

top.sort(
    key=lambda x: float(x["sample_variance"]),
    reverse=True
)

for x in top[:10]:
    print(
        x["orthogroup"],
        "cluster", x["cluster_index"],
        "support=", x["modal_species_support_n"],
        "positive=", x["segment_positive_species_n"],
        "variance=",
        f'{float(x["sample_variance"]):.2f}'
    )


print()
print("Outputs:")
print(OUT_ASSIGN_ALL)
print(OUT_ASSIGN_PRIMARY)
print(OUT_UNITS)
