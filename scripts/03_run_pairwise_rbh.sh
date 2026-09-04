#!/usr/bin/env bash
set -euo pipefail

INDIR="data/refseq"
OUTDIR="work/orthology/rbh"
TMPDIR="work/orthology/tmp"

mkdir -p "$OUTDIR" "$TMPDIR"

FASTAS=()
while IFS= read -r f; do
  FASTAS+=("$f")
done < <(find "$INDIR" -name "protein.faa" | sort)

N=${#FASTAS[@]}

echo "Proteomes found: ${N}"

for ((i=0; i<N; i++)); do
  for ((j=i+1; j<N; j++)); do

    f1="${FASTAS[$i]}"
    f2="${FASTAS[$j]}"

    acc1=$(echo "$f1" | awk -F'/' '{print $3}')
    acc2=$(echo "$f2" | awk -F'/' '{print $3}')

    echo "=== RBH: ${acc1} vs ${acc2} ==="

    mmseqs easy-rbh \
      "$f1" \
      "$f2" \
      "${OUTDIR}/${acc1}__${acc2}.tsv" \
      "${TMPDIR}/${acc1}__${acc2}" \
      --min-seq-id 0.30 \
      -c 0.70 \
      --cov-mode 0 \
      --format-output "query,target,fident,alnlen,qcov,tcov,evalue,bits"

  done
done

echo "=== All pairwise RBH searches finished ==="
