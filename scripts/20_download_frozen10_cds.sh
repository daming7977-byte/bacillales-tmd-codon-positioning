#!/usr/bin/env bash
set -euo pipefail

OUTBASE="data/cds"
mkdir -p "$OUTBASE"

ACCESSIONS=(
  GCF_000009045.1
  GCA_000015785.2
  GCA_000025805.1
  GCF_000011245.1
  GCF_034478925.1
  GCF_006094295.1
  GCF_000017885.4
  GCF_000009785.1
  GCF_034479305.1
  GCF_005671335.1
)

for acc in "${ACCESSIONS[@]}"; do
    echo "========================================"
    echo "Downloading CDS: $acc"
    echo "========================================"

    outdir="$OUTBASE/$acc"
    zipfile="$OUTBASE/${acc}_cds.zip"

    if find "$outdir" -type f -name 'cds_from_genomic.fna' -s 1c 2>/dev/null | grep -q .; then
        echo "$acc already downloaded — skipping."
        continue
    fi

    rm -rf "$outdir"
    rm -f "$zipfile"

    datasets download genome accession "$acc" \
      --include cds \
      --filename "$zipfile"

    mkdir -p "$outdir"
    unzip -q "$zipfile" -d "$outdir"

    cds=$(find "$outdir" -type f -name 'cds_from_genomic.fna' | head -1)

    if [[ -z "$cds" ]]; then
        echo "ERROR: CDS file not found for $acc"
        exit 1
    fi

    echo "CDS file: $cds"
    echo "CDS records: $(grep -c '^>' "$cds")"
done

echo "========================================"
echo "All frozen-panel CDS downloads finished."
echo "========================================"
