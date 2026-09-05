from pathlib import Path
import csv
import shutil

OUT = Path("supplementary/bacillales")
OUT.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# Frozen species panel
# --------------------------------------------------

species = [
    ("GCF_000009045.1", "Bacillus subtilis 168"),
    ("GCA_000015785.2", "Bacillus velezensis FZB42"),
    ("GCA_000025805.1", "Priestia megaterium DSM 319"),
    ("GCF_000011245.1", "Oceanobacillus iheyensis HTE831"),
    ("GCF_034478925.1", "Bacillus licheniformis ATCC 14580"),
    ("GCF_006094295.1", "Bacillus cereus ATCC 14579"),
    ("GCF_000017885.4", "Bacillus pumilus SAFR-032"),
    ("GCF_000009785.1", "Geobacillus kaustophilus HTA426"),
    ("GCF_034479305.1", "Shouchella clausii ATCC 700160"),
    ("GCF_005671335.1", "Halalkalibacterium halodurans LB-1"),
]

with (OUT / "Table_B1_bacillales_species.tsv").open(
    "w", newline=""
) as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["assembly_accession", "organism"])
    w.writerows(species)


# --------------------------------------------------
# Copy frozen analysis tables
# --------------------------------------------------

copies = {
    Path(
        "work/topology/merged/"
        "topology_qualified_families.tsv"
    ):
        OUT / "Table_B2_topology_qualified_families.tsv",

    Path(
        "work/topology/merged/"
        "topology_qualified_tmd_clusters.tsv"
    ):
        OUT / "Table_B3_topology_qualified_tmd_clusters.tsv",

    Path(
        "work/topology/merged/"
        "modal_tmd_primary_set_diagnostic.tsv"
    ):
        OUT / "Table_B4_modal_tmd_primary_set.tsv",

    Path(
        "work/codon/frozen_thresholds/"
        "species_low_adaptation_thresholds.tsv"
    ):
        OUT / "Table_B5_species_low_adaptation_thresholds.tsv",

    Path(
        "work/codon/segments/"
        "low_adaptation_segments.tsv"
    ):
        OUT / "Table_B6_low_adaptation_segments.tsv",

    Path(
        "work/codon/primary_tmd_assignment/"
        "primary_homologous_tmd_unit_summary.tsv"
    ):
        OUT / "Table_B7_primary_tmd_unit_summary.tsv",

    Path(
        "work/codon/permutation/"
        "primary_synonymous_null_1000.tsv"
    ):
        OUT / "Table_B8_synonymous_permutation_null.tsv",

    Path(
        "work/final/"
        "bacillales_replication_final_summary.tsv"
    ):
        OUT / "Table_B9_final_summary.tsv",
}

for src, dst in copies.items():
    if not src.exists():
        raise FileNotFoundError(src)

    shutil.copy2(src, dst)


# --------------------------------------------------
# README
# --------------------------------------------------

readme = """# Bacillales independent analysis

This directory contains compact supplementary tables for the independent
Bacillales analysis accompanying the Enterobacterales primary study.

## Frozen analysis design

- 10 Bacillales genomes
- strict one-to-one orthology
- orthogroups represented in >=8 species
- DeepTMHMM topology prediction
- MAFFT family alignments
- homologous TMD clustering by >=50% overlap of the shorter aligned interval
- cluster support >=80% of family members
- TMD-start MAD <=2 residues
- TMD-end MAD <=2 residues
- no protein contributes multiple TMDs to a retained cluster
- >=3 retained homologous TMD clusters per topology-qualified family
- conserved N-to-C TMD order
- primary proteins required to match the family modal TMD count
- low-adaptation threshold = species-specific bottom decile
- ties retained
- minimum segment length = 3 codons
- initiation and stop codons excluded
- homologous TMD unit support >=8 species
- segment-positive species >=3
- 1,000 within-protein synonymous-codon permutations
- random seed 20260825

## Frozen primary results

TMD start:
- qualifying units: 30
- observed median variance: 2031.10 aa^2
- null median: 1817.125 aa^2
- empirical one-sided P: 0.595405

TMD end:
- qualifying units: 33
- observed median variance: 739.083333 aa^2
- null median: 1693.00 aa^2
- empirical one-sided P: 0.082917

TMD center:
- qualifying units: 32
- observed median variance: 1769.78125 aa^2
- null median: 1813.820833 aa^2
- empirical one-sided P: 0.471528

Interpretation:
The strong Enterobacterales signal was not fully reproduced in Bacillales.
TMD-end-relative variance showed a directional reduction relative to the
synonymous null, whereas TMD-start and TMD-center analyses did not.
"""

(OUT / "README.md").write_text(readme)

print("Created supplementary package:")
print(OUT)

for p in sorted(OUT.iterdir()):
    print(" ", p.name)
