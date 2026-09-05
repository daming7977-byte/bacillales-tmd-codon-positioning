# Bacillales independent analysis

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
