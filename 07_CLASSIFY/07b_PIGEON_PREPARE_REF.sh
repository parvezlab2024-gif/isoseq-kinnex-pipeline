#!/bin/bash
#SBATCH --job-name=pigeon_prep_ref
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 07b: pigeon prepare for the reference (once per annotation).
# Sorts + indexes the GTF (foo.gtf -> foo.sorted.gtf + .pgi), makes the
# genome .fai if missing, and does the same for CAGE_BED if you set one.
# Can run in parallel with collapse.
#
#   cd 07_CLASSIFY
#   sbatch 07b_PIGEON_PREPARE_REF.sh
#
# If this fails, the GTF is almost always the problem. See REFERENCE/README.md.

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 07_CLASSIFY
job_header "pigeon prepare (reference): $BUILD_LABEL"
activate_env
need_tools pigeon

need_file "$ANNOTATION_GTF" "check ANNOTATION_GTF in config.sh"
need_file "$GENOME_FA" "check GENOME_FA in config.sh"
case "$ANNOTATION_GTF" in
    *.sorted.gtf) die "ANNOTATION_GTF should be the unsorted GTF, not the .sorted.gtf" ;;
    *.gtf) ;;
    *) die "ANNOTATION_GTF has to end in .gtf (GFF3 needs converting first, see REFERENCE/README.md)" ;;
esac

# quick sanity checks before pigeon gets to it
head -1 "$ANNOTATION_GTF" | grep -q '^#' && \
    warn "GTF starts with # header lines. GENCODE is fine with this; if prepare fails, strip them (see 07a)."
if ! head -2000 "$ANNOTATION_GTF" | grep -v '^#' | grep -q 'gene_name "'; then
    warn "no gene_name in the first 2000 lines of the GTF. NCBI? run 07a first."
fi
GTF_CHR=$(grep -v '^#' "$ANNOTATION_GTF" | head -1 | cut -f1)
FA_CHR=$(grep -m1 '^>' "$GENOME_FA" | sed 's/^>//; s/[[:space:]].*//')
echo "first chrom in GTF: $GTF_CHR   first chrom in FASTA: $FA_CHR"
if [ -f "${GENOME_FA}.fai" ] && ! cut -f1 "${GENOME_FA}.fai" | grep -qx -- "$GTF_CHR"; then
    die "GTF chromosome '$GTF_CHR' is not in the genome. chr1 vs 1 mismatch? Fix before going on."
fi

echo "pigeon prepare $ANNOTATION_GTF $GENOME_FA"
pigeon prepare "$ANNOTATION_GTF" "$GENOME_FA" --log-level INFO || die "pigeon prepare failed"
need_file "$SORTED_GTF" "pigeon prepare did not make the sorted GTF"
need_file "${SORTED_GTF}.pgi"

if [ -n "$CAGE_BED" ]; then
    need_file "$CAGE_BED" "check CAGE_BED in config.sh"
    first=$(head -1 "$CAGE_BED" | cut -f1)
    case "$first" in
        track*|browser*|\#*|chrom|chromosome)
            die "CAGE_BED has a header line ($first). pigeon wants plain BED, remove it first." ;;
    esac
    cut -f1 "${GENOME_FA}.fai" | grep -qx -- "$first" || \
        die "CAGE_BED chromosome '$first' is not in the genome (chr naming?)"
    echo "pigeon prepare $CAGE_BED"
    pigeon prepare "$CAGE_BED" || die "pigeon prepare failed on CAGE_BED"
    need_file "${CAGE_SORTED}.pgi"
fi

echo ""
echo "ready:"
echo "  $SORTED_GTF"
echo "  ${GENOME_FA}.fai"
[ -n "$CAGE_SORTED" ] && echo "  $CAGE_SORTED"
echo "genes: $(awk -F'\t' '$3=="gene"' "$SORTED_GTF" | wc -l)   transcripts: $(awk -F'\t' '$3=="transcript"' "$SORTED_GTF" | wc -l)"
job_footer
