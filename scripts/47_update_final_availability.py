from pathlib import Path
from shutil import copy2
from docx import Document

SRC = Path(
    "/Users/liming/Desktop/秋山研/2026李明投稿论文/李明博论/"
    "MBE_submission_manuscript_v1.3_submission_candidate.docx"
)

DST = Path(
    "/Users/liming/Desktop/秋山研/2026李明投稿论文/李明博论/"
    "MBE_submission_manuscript_v1.4_final_submission_candidate.docx"
)

copy2(SRC, DST)

doc = Document(DST)

old_data = (
    "Processed analysis tables underlying the main figures and statistical analyses, "
    "including ortholog-group assignments, topology-qualified membrane-protein families, "
    "low-adaptation codon segments, TMD-relative positional measurements, permutation-null "
    "summaries, matched soluble-family controls, and family-level comparison statistics, "
    "are publicly available in the associated reproducibility repository and its archived "
    "Zenodo release (https://github.com/daming7977-byte/enterobacterales-tmd-codon-positioning; "
    "https://doi.org/10.5281/zenodo.22254908). Processed tables underlying the independent "
    "Bacillales analysis are provided as Supplementary Tables B1–B9."
)

new_data = (
    "Processed analysis tables underlying the main figures and statistical analyses, "
    "including ortholog-group assignments, topology-qualified membrane-protein families, "
    "low-adaptation codon segments, TMD-relative positional measurements, permutation-null "
    "summaries, matched soluble-family controls, and family-level comparison statistics, "
    "are publicly available in the Enterobacterales reproducibility repository "
    "(https://github.com/daming7977-byte/enterobacterales-tmd-codon-positioning) and its "
    "archived Zenodo release (https://doi.org/10.5281/zenodo.22254908). Processed tables "
    "underlying the independent Bacillales analysis are provided as Supplementary Tables "
    "B1–B9 and are additionally available in the Bacillales reproducibility repository "
    "(https://github.com/daming7977-byte/bacillales-tmd-codon-positioning), archived at "
    "Zenodo under DOI: 10.5281/zenodo.22345714."
)

old_code = (
    "Custom scripts used for ortholog filtering, topology integration, codon-adaptation "
    "calculation, low-adaptation segment detection, TMD-relative positional analysis, "
    "permutation tests, homologous non-TMD pseudo-anchor controls, matched soluble-family "
    "construction, and family-level statistical analyses are publicly available at "
    "https://github.com/daming7977-byte/enterobacterales-tmd-codon-positioning and are "
    "archived at Zenodo under DOI: 10.5281/zenodo.22254908. Custom scripts for the "
    "independent Bacillales orthology, topology clustering, codon-weight calculation, "
    "low-adaptation segment detection, TMD-relative assignment, and synonymous-permutation "
    "analyses will be deposited with the final reproducibility release."
)

new_code = (
    "Custom scripts used for ortholog filtering, topology integration, codon-adaptation "
    "calculation, low-adaptation segment detection, TMD-relative positional analysis, "
    "permutation tests, homologous non-TMD pseudo-anchor controls, matched soluble-family "
    "construction, and family-level statistical analyses are publicly available at "
    "https://github.com/daming7977-byte/enterobacterales-tmd-codon-positioning and are "
    "archived at Zenodo under DOI: 10.5281/zenodo.22254908. Custom scripts for the "
    "independent Bacillales orthology, topology clustering, codon-weight calculation, "
    "low-adaptation segment detection, TMD-relative assignment, and synonymous-permutation "
    "analyses are publicly available at "
    "https://github.com/daming7977-byte/bacillales-tmd-codon-positioning and are archived "
    "at Zenodo under DOI: 10.5281/zenodo.22345714."
)

changes = 0

for p in doc.paragraphs:
    if p.text == old_data:
        p.text = new_data
        changes += 1
    elif p.text == old_code:
        p.text = new_code
        changes += 1

doc.save(DST)

print("Saved:", DST)
print("Changes:", changes)

if changes != 2:
    print("WARNING: expected exactly 2 paragraph replacements.")
