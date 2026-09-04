#!/usr/bin/env bash
set -euo pipefail

OUTDIR="data/refseq"
mkdir -p "${OUTDIR}"

ACCESSIONS=(
  "GCF_000009045.1"
  "GCA_000025805.1"
  "GCA_000146875.2"
  "GCF_000196035.1"
  "GCF_000013425.1"
  "GCF_000011245.1"
  "GCA_000015785.2"
  "GCF_000007645.1"
)

for ACC in "${ACCESSIONS[@]}"; do
  echo "=== Downloading ${ACC} ==="

  datasets download genome accession "${ACC}" \
    --include protein \
    --filename "${OUTDIR}/${ACC}.zip"

  mkdir -p "${OUTDIR}/${ACC}"

  unzip -o "${OUTDIR}/${ACC}.zip" \
    -d "${OUTDIR}/${ACC}" >/dev/null

done

echo "=== Finished ==="
