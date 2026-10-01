#!/bin/bash
# Step 03: build the FLNC file lists (fofn) used by clustering, collapse and
# per-sample mapping. Quick, run it directly (no sbatch):
#
#   cd 03_MERGE_FOFN
#   bash 03_CREATE_FOFN.sh
#
# Writes to fofn_files/:
#   flnc_ALL.fofn        every refined BAM. This is the one that matters:
#                        cluster once, collapse once, and you get one isoform
#                        set with a count column per sample.
#   flnc_<GROUP>.fofn    one per group in samples.tsv. Optional, for
#                        exploring each group on its own. Do not use these to
#                        build a count matrix for differential analysis;
#                        separate clustering gives isoform sets that are not
#                        comparable to each other.
#   sample_<SAMPLE>.fofn per sample (can be >1 BAM if a sample was on
#                        several cells). Used by 05c.
#   design.tsv           sample -> group, for DTU/DGE tools later.

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "${HERE}/../lib/functions.sh" || exit 1
check_sheet

REFINED="${LIMA_DIR}/results/refined"
mkdir -p "$FOFN_DIR"
rm -f "$FOFN_DIR"/*.fofn

# sort each refined BAM into its sample's list
n=0
for bam in "$REFINED"/*.flnc.bam; do
    [ -s "$bam" ] || continue
    read -r CELL PAIR < <(split_demux_name "$bam")
    s=$(sample_for "$CELL" "$PAIR")
    if [ -z "$s" ]; then
        warn "$(basename "$bam") is not in samples.tsv, leaving it out"
        continue
    fi
    echo "$bam" >> "$FOFN_DIR/sample_${s}.fofn"
    echo "$bam" >> "$FOFN_DIR/flnc_ALL.fofn"
    echo "$bam" >> "$FOFN_DIR/flnc_$(group_of "$s").fofn"
    n=$((n+1))
done
[ "$n" -gt 0 ] || die "no refined BAMs found in $REFINED (run 02b first)"

# every sample in the sheet should have at least one BAM
missing=0
while read -r s; do
    if [ ! -s "$FOFN_DIR/sample_${s}.fofn" ]; then
        warn "no refined BAM for sample $s"
        missing=$((missing+1))
    fi
done < <(sheet_samples)

echo -e "sample\tgroup" > "$FOFN_DIR/design.tsv"
sheet_rows | awk -F'\t' '!s[$2]++{print $2"\t"$3}' >> "$FOFN_DIR/design.tsv"

echo "fofn files in $FOFN_DIR:"
for f in "$FOFN_DIR"/*.fofn; do
    printf "  %-35s %s BAM(s)\n" "$(basename "$f")" "$(wc -l < "$f")"
done
echo ""
cat "$FOFN_DIR/design.tsv"
[ "$missing" -eq 0 ] || die "$missing sample(s) have no data, check the warnings above"
echo ""
echo "next: cd ../04_CLUSTERING && sbatch 04a_RUN_CLUSTER.sh ALL"
