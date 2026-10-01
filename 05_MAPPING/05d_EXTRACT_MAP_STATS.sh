#!/bin/bash
#SBATCH --job-name=map_stats
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 05d (QC): flagstat, MAPQ and per-chromosome counts for every mapped
# BAM that exists (clustered ALL/groups + individual samples).
#
#   cd 05_MAPPING
#   sbatch 05d_EXTRACT_MAP_STATS.sh
#   (when it finishes)  python MAPPING_QC.py
#
# Writes results/qc/<BUILD_LABEL>/mapstats_<DATASET>.tsv

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 05_MAPPING
job_header "mapping stats: $BUILD_LABEL"
activate_env
need_tools samtools

QC="results/qc/${BUILD_LABEL}"
mkdir -p "$QC"

extract() {
    local BAM="$1" DS="$2"
    local OUT="${QC}/mapstats_${DS}.tsv"
    [ -f "$BAM" ] || { echo "  $DS: no BAM, skipping"; return; }
    : > "$OUT"

    # flagstat. Only use lines that look the same in old and new samtools,
    # and work out the primary numbers ourselves.
    local FS TOTAL SEC SUP MAPPED
    FS=$(samtools flagstat -@ 4 "$BAM")
    TOTAL=$(awk '$4=="in" && $5=="total"{print $1; exit}' <<< "$FS")
    SEC=$(awk '$4=="secondary"{print $1; exit}' <<< "$FS")
    SUP=$(awk '$4=="supplementary"{print $1; exit}' <<< "$FS")
    MAPPED=$(awk '$4=="mapped"{print $1; exit}' <<< "$FS")
    TOTAL=${TOTAL:-0}; SEC=${SEC:-0}; SUP=${SUP:-0}; MAPPED=${MAPPED:-0}
    local P=$(( TOTAL - SEC - SUP ))
    local PM=$(( MAPPED - SEC - SUP ))
    {
        printf "flagstat_primary\t%s\n" "$P"
        printf "flagstat_primary_mapped\t%s\n" "$PM"
        printf "flagstat_primary_unmapped\t%s\n" "$(( P - PM ))"
        printf "flagstat_supplementary\t%s\n" "$SUP"
        awk -v a="$PM" -v b="$P" 'BEGIN{ printf "flagstat_map_rate\t%.4f\n", (b>0 ? 100*a/b : 0) }'
    } >> "$OUT"

    # MAPQ of primary alignments
    samtools view -@ 4 -F 0x904 "$BAM" \
        | awk '{c[$5]++} END{for (k in c) printf "mapq_%s\t%d\n", k, c[k]}' | sort -t_ -k2,2n >> "$OUT"

    # reads per chromosome; anything not matching MAIN_CHROM_REGEX -> other
    samtools idxstats "$BAM" \
        | awk -v re="$MAIN_CHROM_REGEX" '$1 != "*" {
              if ($1 ~ re) printf "chr_%s\t%d\n", $1, $3
              else other += $3 }
              END { printf "chr_other\t%d\n", other }' >> "$OUT"

    echo "  $DS: primary=$P mapped=$PM -> $OUT"
}

for G in ALL $(sheet_groups); do
    extract "results/${BUILD_LABEL}/${G}/clustered_${G}.${BUILD_LABEL}.mapped.bam" "$G"
done
for S in $(sheet_samples); do
    extract "results/${BUILD_LABEL}/individual/${S}.${BUILD_LABEL}.mapped.bam" "$S"
done
echo "next: python MAPPING_QC.py"
job_footer
