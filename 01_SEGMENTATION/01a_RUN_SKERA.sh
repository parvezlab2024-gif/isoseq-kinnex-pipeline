#!/bin/bash
#SBATCH --job-name=skera
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 01a: skera split
# Splits the Kinnex/MAS-Seq concatemers (HiFi reads) back into single
# cDNA segments (S-reads). One array task per BAM in HIFI_BAMS.
#
#   cd 01_SEGMENTATION
#   sbatch --array=1-<number of HIFI_BAMS> 01a_RUN_SKERA.sh
#
# Skip this step if your core already gave you segmented reads
# (*.segmented.bam / s-reads). Just drop those in results/segmented/ with
# the name <cell>.segmented.bam and move on to 02.
#
# Check <cell>.segmented.summary.json afterwards. With the right adapter set
# percent_full_array should be well above 90%. For our data mas8 gave 94.5%
# full arrays; the wrong set (mas16) gave 0%.

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 01_SEGMENTATION
job_header "skera split"

activate_env
need_tools skera samtools

INPUT=$(array_pick "${HIFI_BAMS[@]}") || exit 1
CELL=$(cell_id "$INPUT")
mkdir -p results/segmented
OUTPUT="results/segmented/${CELL}.segmented.bam"
LOG="logs/skera.${CELL}.log"

echo "input    : $INPUT"
echo "adapters : $SKERA_ADAPTERS"
echo "output   : $OUTPUT"
need_file "$INPUT" "check HIFI_BAMS in config.sh"
need_file "$SKERA_ADAPTERS" "check SKERA_ADAPTERS in config.sh"

skera split -j "$THREADS_SKERA" "$INPUT" "$SKERA_ADAPTERS" "$OUTPUT" 2> "$LOG" \
    || die "skera failed, see $LOG"
need_file "$OUTPUT"

IN_N=$(count_reads "$INPUT")
OUT_N=$(count_reads "$OUTPUT")
echo ""
echo "HiFi reads in : $IN_N"
echo "S-reads out   : $OUT_N"
awk -v a="$IN_N" -v b="$OUT_N" 'BEGIN{ if (a>0) printf "fold          : %.2fx\n", b/a }'
SUMMARY="${OUTPUT%.bam}.summary.json"
[ -f "$SUMMARY" ] && grep -A2 -E '"(mean_array_size|percent_full_array)"' "$SUMMARY" | grep -E 'id|value'

job_footer
