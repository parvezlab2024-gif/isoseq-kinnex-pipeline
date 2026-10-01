#!/bin/bash
#SBATCH --job-name=pbmm2_index
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 05a: index the genome once (samtools faidx + pbmm2 ISOSEQ index).
# Can run any time before 05b, e.g. while clustering is going.
#
#   cd 05_MAPPING
#   sbatch 05a_INDEX_GENOME.sh
#
# The .mmi is written next to GENOME_FA, so that folder has to be writable.

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 05_MAPPING
job_header "index genome: $BUILD_LABEL"
activate_env
need_tools pbmm2 samtools

case "$GENOME_FA" in
    *.gz) die "GENOME_FA is gzipped. gunzip it first (pbmm2 is fine with .gz but pigeon/faidx want plain or bgzip)." ;;
esac
need_file "$GENOME_FA" "check GENOME_FA in config.sh"

echo "samtools faidx $GENOME_FA"
samtools faidx "$GENOME_FA" || die "faidx failed"
echo "pbmm2 index --preset ISOSEQ -> $GENOME_MMI"
pbmm2 index --preset ISOSEQ -j "$THREADS_MAP" "$GENOME_FA" "$GENOME_MMI" \
    || die "pbmm2 index failed"

echo ""
echo "sequences in genome: $(wc -l < "${GENOME_FA}.fai")"
echo "first few names (these must match the GTF):"
cut -f1 "${GENOME_FA}.fai" | head -5 | sed 's/^/  /'
job_footer
