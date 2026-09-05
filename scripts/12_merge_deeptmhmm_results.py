from pathlib import Path
import json
import csv
from collections import Counter

INDIR = Path("work/topology/deeptmhmm_results")
OUTDIR = Path("work/topology/merged")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_TSV = OUTDIR / "deeptmhmm_all_predictions.tsv"
OUT_SEGMENTS = OUTDIR / "deeptmhmm_all_segments.tsv"

rows = []
segment_rows = []

json_files = sorted(INDIR.glob("batch_*/predictions.json"))
print("Prediction JSON files:", len(json_files))

seen_ids = set()
duplicate_ids = []
skipped_records = []

for jf in json_files:
    batch = jf.parent.name

    with jf.open() as f:
        data = json.load(f)

    print(f"{batch}: raw records = {len(data)}")

    for pred in data:
        seq_id = pred.get("id", "")

        parts = seq_id.split("|", 2)

        # Keep only the sequence records that match our FASTA header format
        if len(parts) != 3:
            skipped_records.append((batch, seq_id))
            continue

        orthogroup, accession, protein_id = parts

        if seq_id in seen_ids:
            duplicate_ids.append(seq_id)
            continue

        seen_ids.add(seq_id)

        protein_type = pred.get("type", "")
        membrane_types = pred.get("membrane_types", [])
        topology_string = pred.get("topology_string", "")
        segments = pred.get("segments", [])

        if isinstance(membrane_types, list):
            membrane_type_text = ";".join(map(str, membrane_types))
        else:
            membrane_type_text = str(membrane_types)

        tm_segments = []

        for seg in segments:
            if not isinstance(seg, list) or len(seg) < 3:
                continue

            label = str(seg[0])
            start = seg[1]
            end = seg[2]

            segment_rows.append({
                "batch": batch,
                "orthogroup": orthogroup,
                "accession": accession,
                "protein_id": protein_id,
                "protein_type": protein_type,
                "segment_label": label,
                "start": start,
                "end": end,
            })

            l = label.lower()

            # DeepTMHMM2 segment labels are expected to include things like:
            # inside, outside, TMhelix, signal, beta, reentrant, interfacial
            # Count only alpha-helical TM segments here.
            is_alpha_tm = (
                l in {"tmhelix", "tm", "transmembrane", "alpha"}
                or "tmhelix" in l
                or "transmembrane" in l
            )

            if is_alpha_tm and "beta" not in l:
                tm_segments.append((start, end))

        rows.append({
            "batch": batch,
            "orthogroup": orthogroup,
            "accession": accession,
            "protein_id": protein_id,
            "protein_type": protein_type,
            "membrane_types": membrane_type_text,
            "n_tmd": len(tm_segments),
            "tmd_coordinates": ";".join(
                f"{s}-{e}" for s, e in tm_segments
            ),
            "topology_string": topology_string,
        })

print()
print("Unique predictions retained:", len(rows))
print("Skipped records:", len(skipped_records))
print("Duplicate IDs:", len(duplicate_ids))
print("Segments loaded:", len(segment_rows))

if skipped_records:
    print("\nFirst skipped records:")
    for x in skipped_records[:20]:
        print(x)

if duplicate_ids:
    print("\nFirst duplicate IDs:")
    for x in duplicate_ids[:20]:
        print(x)

with OUT_TSV.open("w", newline="") as f:
    fieldnames = [
        "batch",
        "orthogroup",
        "accession",
        "protein_id",
        "protein_type",
        "membrane_types",
        "n_tmd",
        "tmd_coordinates",
        "topology_string",
    ]

    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

with OUT_SEGMENTS.open("w", newline="") as f:
    fieldnames = [
        "batch",
        "orthogroup",
        "accession",
        "protein_id",
        "protein_type",
        "segment_label",
        "start",
        "end",
    ]

    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(segment_rows)

print()
print("Output:", OUT_TSV)
print("Segments:", OUT_SEGMENTS)
