#!/bin/bash
# Step 07a (only for NCBI RefSeq GTFs): make the GTF pigeon-friendly.
# GENCODE and Ensembl GTFs do not need this, skip to 07b.
#
#   bash 07a_FIX_NCBI_GTF.sh REFERENCE/GCF_xxx_genomic.gtf
#
# Writes <name>.pigeon.gtf next to the input. Then set ANNOTATION_GTF in
# config.sh to that file.
#
# What it changes (see REFERENCE/README.md for why):
#   - drops # comment lines
#   - NCBI uses gene "xyz" instead of gene_name "xyz" -> adds gene_name
#   - removes the empty transcript_id "" that NCBI puts on gene lines
#   - drops lines that do not have 9 tab separated columns
# We used this on the zebrafish GRCz12tu RefSeq annotation.

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "${HERE}/../lib/functions.sh" || exit 1

IN="${1:-}"
[ -n "$IN" ] || die "usage: bash 07a_FIX_NCBI_GTF.sh <ncbi.gtf>"
IN="$(abspath "$IN")"
need_file "$IN"
case "$IN" in *.gz) die "gunzip the GTF first" ;; esac
OUT="${IN%.gtf}.pigeon.gtf"

echo "in : $IN"
echo "out: $OUT"

# gawk for match() with an array
AWK=$(command -v gawk || echo awk)
"$AWK" 'BEGIN { FS = OFS = "\t" }
    /^#/ { next }
    NF != 9 { next }
    {
        if ($9 !~ /gene_name "/ && match($9, /(^|; )gene "([^"]+)"/, m))
            sub(/gene_id "[^"]*";/, "& gene_name \"" m[2] "\";", $9)
        if ($3 == "gene")
            gsub(/ ?transcript_id "";/, "", $9)
        print
    }' "$IN" > "${OUT}.tmp" || die "awk failed (needs gawk)"
mv "${OUT}.tmp" "$OUT"

echo ""
echo "lines           : $(wc -l < "$OUT")"
echo "with gene_name  : $(grep -c 'gene_name "' "$OUT")"
echo "without         : $(grep -vc 'gene_name "' "$OUT")   (should be 0 or close)"
echo "gene lines with empty transcript_id: $(awk -F'\t' '$3=="gene"' "$OUT" | grep -c 'transcript_id ""')"
echo ""
echo "example:"
awk -F'\t' '$3=="transcript"{print; exit}' "$OUT"
echo ""
echo "now set ANNOTATION_GTF=\"${OUT#${PROJECT_DIR}/}\" in config.sh and run 07b"
