from pathlib import Path
from shutil import copy2
from docx import Document

SRC = Path(
    "/Users/liming/Desktop/秋山研/2026李明投稿论文/李明博论/"
    "MBE_submission_manuscript_v1.4_final_submission_candidate.docx"
)

DST = Path(
    "/Users/liming/Desktop/秋山研/2026李明投稿论文/李明博论/"
    "MBE_submission_manuscript_v1.5_final_submission_candidate.docx"
)

copy2(SRC, DST)

doc = Document(DST)

replacements = [
    # figure
    ("Supplementary Figure B1", "Supplementary Figure S6"),

    # ranges first
    ("Supplementary Tables B1–B9", "Supplementary Tables S7–S15"),
    ("Supplementary Tables B7–B9", "Supplementary Tables S13–S15"),
    ("Supplementary Tables B8–B9", "Supplementary Tables S14–S15"),

    # individual tables
    ("Supplementary Table B9", "Supplementary Table S15"),
    ("Supplementary Table B8", "Supplementary Table S14"),
    ("Supplementary Table B7", "Supplementary Table S13"),
    ("Supplementary Table B6", "Supplementary Table S12"),
    ("Supplementary Table B5", "Supplementary Table S11"),
    ("Supplementary Table B4", "Supplementary Table S10"),
    ("Supplementary Table B3", "Supplementary Table S9"),
    ("Supplementary Table B2", "Supplementary Table S8"),
    ("Supplementary Table B1", "Supplementary Table S7"),

    # plural forms that may occur individually
    ("Supplementary Tables B9", "Supplementary Tables S15"),
    ("Supplementary Tables B8", "Supplementary Tables S14"),
    ("Supplementary Tables B7", "Supplementary Tables S13"),
    ("Supplementary Tables B6", "Supplementary Tables S12"),
    ("Supplementary Tables B5", "Supplementary Tables S11"),
    ("Supplementary Tables B4", "Supplementary Tables S10"),
    ("Supplementary Tables B3", "Supplementary Tables S9"),
    ("Supplementary Tables B2", "Supplementary Tables S8"),
    ("Supplementary Tables B1", "Supplementary Tables S7"),
]

changes = 0

for p in doc.paragraphs:
    original = p.text
    new = original

    for old, repl in replacements:
        new = new.replace(old, repl)

    if new != original:
        p.text = new
        changes += 1

doc.save(DST)

print("Saved:", DST)
print("Paragraphs changed:", changes)
