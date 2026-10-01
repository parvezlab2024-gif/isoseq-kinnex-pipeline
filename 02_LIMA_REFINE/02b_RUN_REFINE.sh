#!/bin/bash
#SBATCH --job-name=refine
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 02b: isoseq refine (--require-polya)
# Trims the poly(A) tail and removes concatemers -> FLNC reads.
# One array task per demuxed BAM (cells x samples). To get the count:
#
#   cd 02_LIMA_REFINE
#   bash 02b_RUN_REFINE.sh --count
#   sbatch --array=1-<that number> 02b_RUN_REFINE.sh
#
# BAMs for barcodes that are not in samples.tsv are skipped (task exits 0).
# This is the output everything downstream uses (results/refined/).

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1

mapfile -t INPUTS < <(demux_bams)
if [ "${1:-}" = "--count" ]; then
    echo "${#INPUTS[@]}"; printf '  %s\n' "${INPUTS[@]}" >&2; exit 0
fi

in_step_dir 02_LIMA_REFINE
job_header "isoseq refine --require-polya"
activate_env
need_tools isoseq samtools
[ "${#INPUTS[@]}" -gt 0 ] || die "no demuxed BAMs in results/demux (run 02a first)"

INPUT=$(array_pick "${INPUTS[@]}") || exit 1
read -r CELL PAIR < <(split_demux_name "$INPUT")
SAMPLE=$(sample_for "$CELL" "$PAIR")
if [ -z "$SAMPLE" ]; then
    echo "$PAIR on $CELL is not in samples.tsv, skipping."
    exit 0
fi

mkdir -p results/refined
OUTPUT="results/refined/${CELL}.${PAIR}.flnc.bam"
LOG="logs/refine.${CELL}.${PAIR}.log"
echo "input  : $INPUT"
echo "sample : $SAMPLE"
echo "output : $OUTPUT"

isoseq refine "$INPUT" "$LIMA_PRIMERS" "$OUTPUT" \
    --require-polya --min-polya-length "$MIN_POLYA_LENGTH" \
    -j "$THREADS_REFINE" --log-level INFO 2> "$LOG" \
    || die "isoseq refine failed, see $LOG"
need_file "$OUTPUT"

# counts come from the refine report, no need to re-read the BAM
REPORT="${OUTPUT%.bam}.filter_summary.report.json"
if [ -f "$REPORT" ]; then
    echo ""
    awk -F': ' '/"name"/ {gsub(/[",]/, "", $2); n = $2}
                /"value"/ {gsub(/[",]/, "", $2); if (n ~ /Reads/) print "  " n ": " $2}' "$REPORT"
fi
job_footer
