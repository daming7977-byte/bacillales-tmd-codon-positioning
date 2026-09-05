#!/usr/bin/env bash
set -euo pipefail

INDIR="work/topology/candidate_fastas"
OUTDIR="work/topology/candidate_alignments"

mkdir -p "$OUTDIR"

for fasta in "$INDIR"/orthogroup_*.faa; do
    name=$(basename "$fasta" .faa)
    outfile="$OUTDIR/${name}.aln.faa"

    echo "Aligning $name"

    if [[ -s "$outfile" ]]; then
        echo "$name already aligned — skipping."
        continue
    fi

    mafft \
      --auto \
      --thread 1 \
      "$fasta" \
      > "$outfile"

done

echo "========================================"
echo "All MAFFT alignments finished."
echo "========================================"
