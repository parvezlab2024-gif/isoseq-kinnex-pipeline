#!/bin/bash
#SBATCH --job-name=lima
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 02a: lima demultiplexing
# Splits each cell's S-reads by cDNA barcode and removes the primers.
# One array task per BAM in HIFI_BAMS.
#
#   cd 02_LIMA_REFINE
#   sbatch --array=1-<number of HIFI_BAMS> 02a_RUN_LIMA.sh
#
# The sample names from samples.tsv are written into the BAM headers here
# (--biosample-csv + --overwrite-biosample-names). isoseq collapse later uses
# those names as the columns of flnc_count.txt. If you skip this, collapse
# puts every read in a single "BioSample_1" column and you lose per-sample
# counts. We learned that one the hard way.
#
# --peek-guess looks at the first reads and drops barcodes that are not
# really there, so unused bc07-bc12 in the primer file are fine.

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 02_LIMA_REFINE
job_header "lima"

activate_env
need_tools lima samtools
check_sheet

HIFI=$(array_pick "${HIFI_BAMS[@]}") || exit 1
CELL=$(cell_id "$HIFI")
INPUT="${SEG_DIR}/results/segmented/${CELL}.segmented.bam"
mkdir -p results/demux
PREFIX="results/demux/${CELL}"
LOG="logs/lima.${CELL}.log"
need_file "$INPUT" "run 01a first"
need_file "$LIMA_PRIMERS" "check LIMA_PRIMERS in config.sh"

# biosample CSV for this cell only
BIOCSV="results/demux/${CELL}.biosamples.csv"
echo "Barcodes,Bio Sample" > "$BIOCSV"
sheet_rows | awk -F'\t' -v cell="$CELL" '$4=="all" || index(cell,$4)>0 {print $1","$2}' >> "$BIOCSV"
[ "$(wc -l < "$BIOCSV")" -gt 1 ] || die "no samples.tsv rows apply to cell $CELL"

echo "input    : $INPUT"
echo "primers  : $LIMA_PRIMERS"
echo "samples  :"; sed 's/^/    /' "$BIOCSV"

lima "$INPUT" "$LIMA_PRIMERS" "${PREFIX}.bam" \
    --isoseq --peek-guess \
    --biosample-csv "$BIOCSV" --overwrite-biosample-names \
    -j "$THREADS_LIMA" --log-level INFO 2> "$LOG" \
    || die "lima failed, see $LOG"

echo ""
echo "demuxed BAMs for $CELL:"
for bam in "${PREFIX}".*.bam; do
    [ -f "$bam" ] || continue
    pair=$(split_demux_name "$bam" | cut -f2)
    s=$(sample_for "$CELL" "$pair")
    printf "  %-40s %-15s %s reads\n" "$pair" "${s:-NOT_IN_SHEET}" "$(count_reads "$bam")"
done
echo "summary: ${PREFIX}.lima.summary"
job_footer
