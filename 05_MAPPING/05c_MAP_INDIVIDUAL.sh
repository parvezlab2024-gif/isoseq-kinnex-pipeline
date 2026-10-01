#!/bin/bash
#SBATCH --job-name=map_sample
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 05c (optional): map each sample's FLNC reads (not clustered) to the
# genome. These are for IGV / looking at a locus sample by sample. Collapse
# does not use them, but the mapping QC will include them if present.
# One array task per sample in samples.tsv:
#
#   cd 05_MAPPING
#   bash 05c_MAP_INDIVIDUAL.sh --count
#   sbatch --array=1-<that number> 05c_MAP_INDIVIDUAL.sh

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1

mapfile -t SAMPLES < <(sheet_samples)
if [ "${1:-}" = "--count" ]; then
    echo "${#SAMPLES[@]}"; printf '  %s\n' "${SAMPLES[@]}" >&2; exit 0
fi

in_step_dir 05_MAPPING
SAMPLE=$(array_pick "${SAMPLES[@]}") || exit 1
job_header "pbmm2 align (sample $SAMPLE) -> $BUILD_LABEL"
activate_env
need_tools pbmm2 samtools

FOFN="${FOFN_DIR}/sample_${SAMPLE}.fofn"
need_file "$FOFN" "run 03_CREATE_FOFN.sh first"
need_file "$GENOME_MMI" "run 05a first"

OUTDIR="results/${BUILD_LABEL}/individual"
mkdir -p "$OUTDIR"
OUTPUT="${OUTDIR}/${SAMPLE}.${BUILD_LABEL}.mapped.bam"
LOG="${OUTDIR}/${SAMPLE}.pbmm2.log"
echo "input : $FOFN ($(wc -l < "$FOFN") BAMs)"
echo "output: $OUTPUT"

pbmm2 align --preset ISOSEQ --sort --unmapped -j "$THREADS_MAP" --log-level INFO \
    "$GENOME_MMI" "$FOFN" "$OUTPUT" 2> "$LOG" \
    || die "pbmm2 failed, see $LOG"
need_file "$OUTPUT"
samtools index "$OUTPUT"
echo "primary mapped: $(samtools view -c -F 0x904 "$OUTPUT")   unmapped: $(samtools view -c -f 4 "$OUTPUT")"
job_footer
