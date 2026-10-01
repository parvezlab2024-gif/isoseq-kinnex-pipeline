#!/bin/bash
#SBATCH --job-name=collapse
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 06a: isoseq collapse
# Merges mapped transcripts into unique isoforms and counts FLNC reads per
# isoform per sample (flnc_count.txt). For ALL, flnc_count.txt is your
# isoform x sample count matrix.
#
#   cd 06_COLLAPSE
#   sbatch 06a_COLLAPSE.sh ALL
#   sbatch 06a_COLLAPSE.sh <GROUP>     (only if you ran 04a/05b for it)
#
# COLLAPSE_OPTS is in config.sh. We used --do-not-collapse-extra-5exons so
# 5' truncated models are kept apart from the full-length one.

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 06_COLLAPSE

GROUP="${1:-ALL}"
check_group "$GROUP"
job_header "isoseq collapse: $GROUP ($BUILD_LABEL)"
activate_env
need_tools isoseq

MAPPED="${MAP_DIR}/results/${BUILD_LABEL}/${GROUP}/clustered_${GROUP}.${BUILD_LABEL}.mapped.bam"
FOFN="${FOFN_DIR}/flnc_${GROUP}.fofn"
need_file "$MAPPED" "run 05b $GROUP first"
need_file "$FOFN" "run 03 first"

OUTDIR="results/${BUILD_LABEL}/${GROUP}"
mkdir -p "$OUTDIR"
OUT_GFF="${OUTDIR}/clustered_${GROUP}.${BUILD_LABEL}.collapsed.gff"
LOG="${OUTDIR}/collapse.log"
echo "mapped : $MAPPED"
echo "flnc   : $FOFN"
echo "output : $OUT_GFF"

isoseq collapse ${COLLAPSE_OPTS} -j "$THREADS_COLLAPSE" --log-level INFO \
    "$MAPPED" "$FOFN" "$OUT_GFF" 2> "$LOG" \
    || die "collapse failed, see $LOG"
need_file "$OUT_GFF"

COUNTS="${OUT_GFF%.gff}.flnc_count.txt"
echo "isoforms: $(grep -c $'\ttranscript\t' "$OUT_GFF")"
if [ -s "$COUNTS" ]; then
    echo "flnc_count.txt columns: $(head -1 "$COUNTS")"
    # a single BioSample_1 column means the lima biosample names were lost
    head -1 "$COUNTS" | grep -q "BioSample_1" && \
        warn "flnc_count.txt has a BioSample_1 column. Sample names did not make it through lima (see 02a)."
else
    warn "no flnc_count.txt, per-sample counts are missing"
fi
job_footer
