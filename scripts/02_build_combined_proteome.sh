#!/usr/bin/env bash
set -euo pipefail

OUT="work/orthology/combined_8species.faa"
: > "$OUT"

for f in $(find data/refseq -name "protein.faa" | sort); do
    acc=$(echo "$f" | awk -F'/' '{print $3}')

    awk -v acc="$acc" '
    /^>/ {
        sub(/^>/, "")
        print ">" acc "|" $0
        next
    }
    {print}
    ' "$f" >> "$OUT"
done

echo "Combined FASTA: $OUT"
echo -n "Total proteins: "
grep -c "^>" "$OUT"
