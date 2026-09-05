from docx import Document

DOCX = "/Users/liming/Desktop/秋山研/2026李明投稿论文/李明博论/MBE_submission_manuscript_v1.2.docx"

doc = Document(DOCX)

fixes = {
    "153topology-qualified": "153 topology-qualified",
    "TMD countof": "TMD count of",
    "speciesafter": "species after",
    "1,372proteins": "1,372 proteins",
    "Low-adaptationsegments": "Low-adaptation segments",
    "the segmentwas": "the segment was",
    "species, anda unit": "species, and a unit",
    "pseudo-anchorcontrols": "pseudo-anchor controls",
    "publicly availableat": "publicly available at",
}

changes = 0

for p in doc.paragraphs:
    old = p.text
    new = old

    for a, b in fixes.items():
        new = new.replace(a, b)

    if new != old:
        p.text = new
        changes += 1

doc.save(DOCX)

print("Saved:", DOCX)
print("Paragraphs corrected:", changes)
