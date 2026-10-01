#!/bin/bash
#SBATCH --job-name=map_clustered
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 05b: map clustered transcripts to the genome (feeds collapse).
#
#   cd 05_MAPPING
#   sbatch 05b_MAP_CLUSTERED.sh ALL
#   sbatch 05b_MAP_CLUSTERED.sh <GROUP>     (only if you clustered that group)
#
# --unmapped keeps unmapped records in the BAM so the QC can count them.

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 05_MAPPING

GROUP="${1:-ALL}"
check_group "$GROUP"
job_header "pbmm2 align (clustered $GROUP) -> $BUILD_LABEL"
activate_env
need_tools pbmm2 samtools

INPUT="${CLUST_DIR}/results/${GROUP}/clustered_${GROUP}.bam"
need_file "$INPUT" "run 04a for $GROUP first"
need_file "$GENOME_MMI" "run 05a first"

OUTDIR="results/${BUILD_LABEL}/${GROUP}"
mkdir -p "$OUTDIR"
OUTPUT="${OUTDIR}/clustered_${GROUP}.${BUILD_LABEL}.mapped.bam"
LOG="${OUTDIR}/pbmm2.log"
echo "input : $INPUT"
echo "index : $GENOME_MMI"
echo "output: $OUTPUT"

pbmm2 align --preset ISOSEQ --sort --unmapped -j "$THREADS_MAP" --log-level INFO \
    "$GENOME_MMI" "$INPUT" "$OUTPUT" 2> "$LOG" \
    || die "pbmm2 failed, see $LOG"
need_file "$OUTPUT"
samtools index "$OUTPUT"

IN_N=$(count_reads "$INPUT")
MAPPED=$(samtools view -c -F 0x904 "$OUTPUT")
echo "transcripts in: $IN_N   primary mapped: $MAPPED"
awk -v a="$MAPPED" -v b="$IN_N" 'BEGIN{ if (b>0) printf "mapping rate: %.2f%%\n", 100*a/b }'
job_footer
