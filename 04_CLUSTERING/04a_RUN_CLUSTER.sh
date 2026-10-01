#!/bin/bash
#SBATCH --job-name=cluster2
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 04a: isoseq cluster2
#
#   cd 04_CLUSTERING
#   sbatch 04a_RUN_CLUSTER.sh ALL          <- the one you need
#   sbatch 04a_RUN_CLUSTER.sh <GROUP>      <- optional, one group only
#
# Memory: pooling 6 samples / ~73 M FLNC reads needed 256G for us.
# One group of 3 fit in 128G. If the job dies with exit 137 it ran out of
# memory; bump --mem and resubmit.

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 04_CLUSTERING

GROUP="${1:-ALL}"
check_group "$GROUP"
job_header "isoseq cluster2: $GROUP"
activate_env
need_tools isoseq samtools

FOFN="${FOFN_DIR}/flnc_${GROUP}.fofn"
need_file "$FOFN" "run 03_CREATE_FOFN.sh first"
while read -r f; do need_file "$f"; done < "$FOFN"
echo "input: $FOFN ($(wc -l < "$FOFN") BAMs)"

OUTDIR="results/${GROUP}"
mkdir -p "$OUTDIR"
OUTPUT="${OUTDIR}/clustered_${GROUP}.bam"
LOG="logs/cluster.${GROUP}.log"

isoseq cluster2 "$FOFN" "$OUTPUT" ${CLUSTER_OPTS} \
    -j "$THREADS_CLUSTER" --log-level INFO 2> "$LOG"
rc=$?
if [ $rc -ne 0 ]; then
    [ $rc -eq 137 ] && echo "exit 137 = out of memory, raise --mem"
    die "cluster2 failed (exit $rc), see $LOG"
fi
need_file "$OUTPUT"

echo "clustered transcripts: $(count_reads "$OUTPUT")"
echo "output: $OUTPUT"
job_footer
