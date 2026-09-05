from pathlib import Path
import csv
import statistics
import hashlib

BASE = Path(".")
OUTDIR = Path("work/final")
OUTDIR.mkdir(parents=True, exist_ok=True)

PERM = Path(
    "work/codon/permutation/"
    "primary_synonymous_null_1000.tsv"
)

UNITS = Path(
    "work/codon/primary_tmd_assignment/"
    "primary_homologous_tmd_unit_summary.tsv"
)

THRESHOLDS = Path(
    "work/codon/frozen_thresholds/"
    "species_low_adaptation_thresholds.tsv"
)

SEGMENTS = Path(
    "work/codon/segments/"
    "low_adaptation_segments.tsv"
)

SUMMARY_TSV = OUTDIR / "bacillales_replication_final_summary.tsv"
SUMMARY_MD = OUTDIR / "bacillales_replication_final_summary.md"
HASHES = OUTDIR / "bacillales_replication_file_hashes.tsv"


# --------------------------------------------------
# observed statistics
# --------------------------------------------------

observed = {}

with UNITS.open() as f:
    r = csv.DictReader(f, delimiter="\t")

    by_anchor = {
        "start": [],
        "end": [],
        "center": [],
    }

    for x in r:
        if x["qualifying_unit"] != "1":
            continue

        by_anchor[x["anchor_type"]].append(
            float(x["sample_variance"])
        )

for anchor, vals in by_anchor.items():
    observed[anchor] = {
        "qualifying_units": len(vals),
        "median_variance": statistics.median(vals),
    }


# --------------------------------------------------
# permutation null
# --------------------------------------------------

with PERM.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

final = {}

for anchor in ["start", "end", "center"]:

    vals = [
        float(x[f"{anchor}_median_variance"])
        for x in rows
    ]

    units = [
        int(x[f"{anchor}_qualifying_units"])
        for x in rows
    ]

    obs = observed[anchor]["median_variance"]

    k = sum(x <= obs for x in vals)

    p = (k + 1) / (len(vals) + 1)

    final[anchor] = {
        "observed_units":
            observed[anchor]["qualifying_units"],

        "observed_median_variance":
            obs,

        "null_units_median":
            statistics.median(units),

        "null_median_variance":
            statistics.median(vals),

        "null_le_observed":
            k,

        "permutations":
            len(vals),

        "empirical_p":
            p,
    }


# --------------------------------------------------
# basic dataset counts
# --------------------------------------------------

segment_n = 0
segment_proteins = set()
segment_families = set()

with SEGMENTS.open() as f:
    r = csv.DictReader(f, delimiter="\t")

    for x in r:
        segment_n += 1
        segment_proteins.add(x["sequence_id"])
        segment_families.add(x["orthogroup"])


threshold_species = 0

with THRESHOLDS.open() as f:
    r = csv.DictReader(f, delimiter="\t")
    threshold_species = sum(1 for _ in r)


# --------------------------------------------------
# write TSV
# --------------------------------------------------

with SUMMARY_TSV.open("w", newline="") as f:

    fields = [
        "anchor",
        "observed_units",
        "observed_median_variance",
        "null_units_median",
        "null_median_variance",
        "null_le_observed",
        "permutations",
        "empirical_p",
    ]

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()

    for anchor in ["start", "end", "center"]:
        row = {"anchor": anchor}
        row.update(final[anchor])
        w.writerow(row)


# --------------------------------------------------
# write markdown record
# --------------------------------------------------

with SUMMARY_MD.open("w") as f:

    f.write("# Bacillales replication final frozen summary\n\n")

    f.write("## Dataset\n\n")
    f.write("- Species: 10\n")
    f.write("- Topology-qualified families before modal filter: 153\n")
    f.write("- Qualified-family proteins: 1432\n")
    f.write("- Modal-TMD-matching proteins: 1372\n")
    f.write("- Frozen qualified homologous TMD clusters: 1086\n")
    f.write("- Primary clusters represented in >=8 modal-matching species: 994\n")
    f.write(f"- Low-adaptation segments: {segment_n}\n")
    f.write(f"- Proteins with >=1 low-adaptation segment: {len(segment_proteins)}\n")
    f.write(f"- Families with >=1 low-adaptation segment: {len(segment_families)}\n")
    f.write(f"- Species-specific frozen thresholds: {threshold_species}\n\n")

    f.write("## Frozen primary parameters\n\n")
    f.write("- Low-adaptation threshold: species-specific bottom decile\n")
    f.write("- Ties at threshold retained\n")
    f.write("- Minimum low-adaptation run: 3 consecutive codons\n")
    f.write("- Initiation codon excluded\n")
    f.write("- Stop codons excluded\n")
    f.write("- Primary proteins: observed TMD count equals family modal TMD count\n")
    f.write("- Homologous TMD cluster support: >=8 species after modal filter\n")
    f.write("- Qualifying unit: >=3 segment-positive species\n")
    f.write("- Synonymous permutation: within protein and amino-acid class\n")
    f.write("- Permutations: 1000\n")
    f.write("- Random seed: 20260825\n\n")

    f.write("## Final observed versus synonymous-null results\n\n")
    f.write(
        "| Anchor | Qualifying units | Observed median variance | "
        "Null median variance | Empirical P |\n"
    )
    f.write("|---|---:|---:|---:|---:|\n")

    for anchor in ["start", "end", "center"]:
        x = final[anchor]

        f.write(
            f"| {anchor} | "
            f"{x['observed_units']} | "
            f"{x['observed_median_variance']:.6f} | "
            f"{x['null_median_variance']:.6f} | "
            f"{x['empirical_p']:.6f} |\n"
        )

    f.write("\n## Interpretation frozen before manuscript revision\n\n")
    f.write(
        "- TMD-start-relative positional constraint did not replicate "
        "in Bacillales.\n"
    )
    f.write(
        "- TMD-end-relative variance was directionally lower than the "
        "synonymous null, but did not reach permutation significance.\n"
    )
    f.write(
        "- TMD-center analysis did not show evidence of positional constraint.\n"
    )
    f.write(
        "- The Bacillales analysis is therefore interpreted as partial, "
        "anchor-specific support rather than full replication.\n"
    )
    f.write(
        "- No primary thresholds, topology criteria, segment definitions, "
        "species panel, or assignment rules will be modified after inspection "
        "of these results.\n"
    )


# --------------------------------------------------
# hashes of key final files
# --------------------------------------------------

files = [
    THRESHOLDS,
    SEGMENTS,
    UNITS,
    PERM,
    SUMMARY_TSV,
    SUMMARY_MD,
]

with HASHES.open("w", newline="") as f:

    w = csv.writer(f, delimiter="\t")
    w.writerow(["file", "sha256"])

    for path in files:

        h = hashlib.sha256()

        with path.open("rb") as fh:
            for block in iter(
                lambda: fh.read(1024 * 1024),
                b""
            ):
                h.update(block)

        w.writerow([
            str(path),
            h.hexdigest()
        ])


print("Frozen summary:", SUMMARY_TSV)
print("Frozen markdown:", SUMMARY_MD)
print("Hashes:", HASHES)

print()
for anchor in ["start", "end", "center"]:
    x = final[anchor]

    print(
        anchor,
        "units=", x["observed_units"],
        "observed=", f'{x["observed_median_variance"]:.6f}',
        "null=", f'{x["null_median_variance"]:.6f}',
        "P=", f'{x["empirical_p"]:.6f}',
    )
