from docx import Document

DOCX = "/Users/liming/Desktop/秋山研/2026李明投稿论文/李明博论/MBE_submission_manuscript_v1.2.docx"

doc = Document(DOCX)

replacements = {
    "analysisacross": "analysis across",
    "weightscalculated": "weights calculated",
    "examinedas": "examined as",
}

changes = 0

for p in doc.paragraphs:
    old = p.text
    new = old

    for a, b in replacements.items():
        new = new.replace(a, b)

    if new != old:
        p.text = new
        changes += 1

doc.save(DOCX)

print("Paragraphs corrected:", changes)
print("Saved:", DOCX)
