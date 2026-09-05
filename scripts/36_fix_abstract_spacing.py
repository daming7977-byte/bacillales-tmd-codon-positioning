from docx import Document

DOCX = "/Users/liming/Desktop/秋山研/2026李明投稿论文/李明博论/MBE_submission_manuscript_v1.2.docx"

doc = Document(DOCX)

replacements = {
    "unclear.Here": "unclear. Here",
    "quantified asthe": "quantified as the",
    "first TMDs,under": "first TMDs, under",
}

n = 0

for p in doc.paragraphs:
    old = p.text
    new = old

    for a, b in replacements.items():
        new = new.replace(a, b)

    if new != old:
        p.text = new
        n += 1

doc.save(DOCX)

print("Paragraphs corrected:", n)
print("Saved:", DOCX)
