#!/bin/bash
#SBATCH --job-name=pigeon_classify
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 07c: pigeon prepare (collapsed GFF) -> classify -> filter
#
#   cd 07_CLASSIFY
#   sbatch 07c_PIGEON_CLASSIFY.sh ALL
#   sbatch 07c_PIGEON_CLASSIFY.sh <GROUP>     (optional)
#
# Needs 06a for this group and 07b done.
# If CAGE_BED / POLYA_MOTIFS are set in config.sh they get passed to classify
# (fills the within_CAGE_peak / polyA_motif columns).
#
# Main outputs in results/<BUILD_LABEL>/<GROUP>/:
#   clustered_<GROUP>_classification.txt              everything
#   clustered_<GROUP>_classification.filtered_lite_classification.txt
#   clustered_<GROUP>.sorted.filtered_lite.gff        isoforms that passed
#   *.summary.txt                                     human-readable counts

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 07_CLASSIFY

GROUP="${1:-ALL}"
check_group "$GROUP"
job_header "pigeon classify + filter: $GROUP ($BUILD_LABEL)"
activate_env
need_tools pigeon

SRC="${COLL_DIR}/results/${BUILD_LABEL}/${GROUP}"
COLLAPSED="${SRC}/clustered_${GROUP}.${BUILD_LABEL}.collapsed.gff"
FLNC_COUNT="${SRC}/clustered_${GROUP}.${BUILD_LABEL}.collapsed.flnc_count.txt"
need_file "$COLLAPSED" "run 06a $GROUP first"
need_file "$FLNC_COUNT" "run 06a $GROUP first"
need_file "$SORTED_GTF" "run 07b first"
need_file "${SORTED_GTF}.pgi" "run 07b first"
need_file "${GENOME_FA}.fai" "run 05a or 07b first"

OUTDIR="results/${BUILD_LABEL}/${GROUP}"
mkdir -p "$OUTDIR"
PREFIX="clustered_${GROUP}"

# 1. sort/index the collapsed GFF. Work on a copy so 06_COLLAPSE stays as is.
cp "$COLLAPSED" "${OUTDIR}/${PREFIX}.gff"
pigeon prepare "${OUTDIR}/${PREFIX}.gff" || die "pigeon prepare failed on the collapsed GFF"
SORTED_GFF="${OUTDIR}/${PREFIX}.sorted.gff"
need_file "$SORTED_GFF"
rm -f "${OUTDIR}/${PREFIX}.gff"

# 2. classify
EXTRA=()
if [ -n "$CAGE_BED" ]; then
    need_file "${CAGE_SORTED}.pgi" "CAGE_BED is set but not prepared, rerun 07b"
    EXTRA+=(--cage-peak "$CAGE_SORTED")
fi
if [ -n "$POLYA_MOTIFS" ]; then
    need_file "$POLYA_MOTIFS"
    EXTRA+=(--poly-a "$POLYA_MOTIFS")
fi
echo "pigeon classify ${SORTED_GFF} ${SORTED_GTF} ${GENOME_FA} --flnc ${FLNC_COUNT} ${EXTRA[*]:-}"
pigeon classify "$SORTED_GFF" "$SORTED_GTF" "$GENOME_FA" \
    --flnc "$FLNC_COUNT" ${EXTRA[@]+"${EXTRA[@]}"} \
    -d "$OUTDIR" -o "$PREFIX" -j "$THREADS_PIGEON" --log-level INFO \
    || die "pigeon classify failed"
CLASS="${OUTDIR}/${PREFIX}_classification.txt"
need_file "$CLASS"

# 3. filter (intra-priming, RT switching, low coverage/non-canonical)
pigeon filter "$CLASS" --isoforms "$SORTED_GFF" -j "$THREADS_PIGEON" --log-level INFO \
    || die "pigeon filter failed"

echo ""
echo "--- classify summary ---"
sed -n '1,/^Junctions/p' "${OUTDIR}/${PREFIX}.summary.txt" 2>/dev/null
echo "--- after filter ---"
sed -n '1,4p' "${OUTDIR}/${PREFIX}_classification.filtered.summary.txt" 2>/dev/null
echo ""
ls -1 "$OUTDIR" | sed 's/^/  /'
job_footer
