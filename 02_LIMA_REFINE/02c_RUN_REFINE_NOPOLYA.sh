#!/bin/bash
#SBATCH --job-name=refine_nopolya
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 02c (optional, QC only): isoseq refine WITHOUT --require-polya
# Same as 02b but keeps reads with no poly(A) tail. Only used by REFINE_QC.py
# to show how many reads the poly(A) requirement throws away. A big gap
# points to internal priming or 3' degradation. Nothing downstream uses it.
#
#   cd 02_LIMA_REFINE
#   sbatch --array=1-$(bash 02b_RUN_REFINE.sh --count 2>/dev/null) 02c_RUN_REFINE_NOPOLYA.sh

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 02_LIMA_REFINE
job_header "isoseq refine (no poly(A) requirement)"
activate_env
need_tools isoseq

mapfile -t INPUTS < <(demux_bams)
[ "${#INPUTS[@]}" -gt 0 ] || die "no demuxed BAMs in results/demux (run 02a first)"
INPUT=$(array_pick "${INPUTS[@]}") || exit 1
read -r CELL PAIR < <(split_demux_name "$INPUT")
if [ -z "$(sample_for "$CELL" "$PAIR")" ]; then
    echo "$PAIR on $CELL is not in samples.tsv, skipping."
    exit 0
fi

mkdir -p results/refined_nopolya
OUTPUT="results/refined_nopolya/${CELL}.${PAIR}.flnc.bam"
LOG="logs/refine_nopolya.${CELL}.${PAIR}.log"

isoseq refine "$INPUT" "$LIMA_PRIMERS" "$OUTPUT" \
    --min-polya-length "$MIN_POLYA_LENGTH" \
    -j "$THREADS_REFINE" --log-level INFO 2> "$LOG" \
    || die "isoseq refine failed, see $LOG"
echo "wrote $OUTPUT"
job_footer
