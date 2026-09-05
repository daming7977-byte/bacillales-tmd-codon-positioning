from docx import Document

DOCX = "/Users/liming/Desktop/秋山研/2026李明投稿论文/李明博论/MBE_submission_manuscript_v1.2.docx"

doc = Document(DOCX)

terms = [
    "Bacillales",
    "Enterobacterales",
    "replicat",
    "generaliz",
    "univers",
    "phylogenetic",
    "0.0829",
    "739.08",
    "1,693",
]

for term in terms:
    print()
    print("=" * 80)
    print("SEARCH:", term)
    print("=" * 80)

    n = 0

    for i, p in enumerate(doc.paragraphs):
        if term.lower() in p.text.lower():
            n += 1
            print(f"[Paragraph {i}]")
            print(p.text)
            print()

    print("Matches:", n)
