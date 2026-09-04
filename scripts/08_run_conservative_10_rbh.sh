#!/usr/bin/env bash
set -euo pipefail

PANEL="config/panel_conservative_10.tsv"
OUTDIR="work/orthology/rbh_conservative10"
TMPDIR="work/orthology/tmp_conservative10"

mkdir -p "$OUTDIR" "$TMPDIR"

ACCESSIONS=()

while IFS=$'\t' read -r accession species strain; do
    [[ "$accession" == "accession" ]] && continue
    ACCESSIONS+=("$accession")
done < "$PANEL"

N=${#ACCESSIONS[@]}
echo "Panel species: $N"

for ((i=0; i<N; i++)); do
  for ((j=i+1; j<N; j++)); do

    acc1="${ACCESSIONS[$i]}"
    acc2="${ACCESSIONS[$j]}"

    f1="data/refseq/${acc1}/ncbi_dataset/data/${acc1}/protein.faa"
    f2="data/refseq/${acc2}/ncbi_dataset/data/${acc2}/protein.faa"

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

echo "=== Conservative 10-species RBH finished ==="
