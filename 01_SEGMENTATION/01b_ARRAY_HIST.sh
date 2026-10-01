#!/bin/bash
#SBATCH --job-name=array_hist
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 01b (optional, QC only): exact segments-per-HiFi-read histogram
# Needs one full pass over the segmented BAM (ours was ~65 GB), so run it as
# a job, not on the login node. SKERA_QC.py uses the table if it exists,
# otherwise it estimates the distribution from a subsample.
#
#   cd 01_SEGMENTATION
#   sbatch --array=1-<number of HIFI_BAMS> 01b_ARRAY_HIST.sh

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 01_SEGMENTATION
job_header "segments per read histogram"

activate_env
need_tools samtools

HIFI=$(array_pick "${HIFI_BAMS[@]}") || exit 1
CELL=$(cell_id "$HIFI")
BAM="results/segmented/${CELL}.segmented.bam"
OUT="results/qc/${CELL}_array_hist.tsv"
mkdir -p results/qc
need_file "$BAM" "run 01a first"

# S-reads from the same ZMW (2nd field of the read name) come out of skera
# next to each other, so a run-length count gives segments per HiFi read.
echo -e "segments\tnum_reads" > "$OUT"
samtools view -@ 4 "$BAM" \
  | awk -F'\t' '{
        split($1, a, "/"); z = a[2]
        if (NR > 1 && z != p) { h[c]++; c = 0 }
        c++; p = z
    } END { if (c > 0) h[c]++; for (k in h) print k"\t"h[k] }' \
  | sort -n >> "$OUT"

cat "$OUT"
echo "reads   : $(awk 'NR>1{s+=$2} END{print s}' "$OUT")  (should match 'reads' in the summary.json)"
echo "S-reads : $(awk 'NR>1{s+=$1*$2} END{print s}' "$OUT")  (should match 's_reads')"
job_footer
