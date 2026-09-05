from docx import Document

DOCX = "/Users/liming/Desktop/秋山研/2026李明投稿论文/李明博论/MBE_submission_manuscript_v1.2.docx"

OLD = (
    "These results reveal evolutionarily constrained positioning of "
    "low-adaptation synonymous codon segments relative to homologous TMD "
    "boundaries in Enterobacterales, while also defining the phylogenetic "
    "limits of this pattern and leaving open whether it directly reflects "
    "modulation of translation kinetics, membrane insertion, cotranslational "
    "folding, or another selective constraint."
)

NEW = (
    "These results reveal evolutionarily constrained positioning of "
    "low-adaptation synonymous codon segments relative to homologous TMD "
    "boundaries in Enterobacterales, while the Bacillales analysis suggests "
    "that the strength and anchor dependence of this pattern may vary among "
    "bacterial lineages. Whether the observed positional constraint directly "
    "reflects modulation of translation kinetics, membrane insertion, "
    "cotranslational folding, or another selective pressure remains unresolved."
)

doc = Document(DOCX)

found = False

for p in doc.paragraphs:
    if OLD in p.text:
        p.text = p.text.replace(OLD, NEW)
        found = True
        break

if not found:
    raise RuntimeError("Expected Introduction sentence not found.")

doc.save(DOCX)

print("Introduction wording refined.")
