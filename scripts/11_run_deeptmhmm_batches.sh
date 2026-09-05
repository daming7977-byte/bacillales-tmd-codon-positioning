#!/usr/bin/env bash
set -euo pipefail

INDIR="work/topology/deeptmhmm_batches"
OUTBASE="work/topology/deeptmhmm_results"

mkdir -p "$OUTBASE"

for fasta in "$INDIR"/batch_*.faa; do
    name=$(basename "$fasta" .faa)
    outdir="$OUTBASE/$name"

    echo "========================================"
    echo "Running $name"
    echo "Input:  $fasta"
    echo "Output: $outdir"
    echo "========================================"

    # 已完成的 batch 自动跳过
    if [[ -s "$outdir/predictions.json" ]]; then
        echo "$name already completed — skipping."
        continue
    fi

    rm -rf "$outdir"

    dtm2 \
      "$fasta" \
      "$outdir" \
      --device cpu \
      --batch-size 4

    echo "$name finished."
done

echo "========================================"
echo "All DeepTMHMM batches finished."
echo "========================================"
