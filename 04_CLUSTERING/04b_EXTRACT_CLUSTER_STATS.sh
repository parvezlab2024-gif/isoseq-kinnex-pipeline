#!/bin/bash
#SBATCH --job-name=cluster_stats
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 04b (QC): pull numbers out of the cluster2 output for CLUSTER_QC.py.
# Runs over ALL plus any group that was clustered.
#
#   cd 04_CLUSTERING
#   sbatch 04b_EXTRACT_CLUSTER_STATS.sh
#   (when it finishes)  python CLUSTER_QC.py
#
# Writes results/qc/cluster_stats_<GROUP>.tsv (key <tab> value):
#   flnc_input            FLNC reads going in (from the refine reports)
#   total_isoforms, support_*   reads per cluster, from cluster_report.csv
#   len_*                 transcript lengths, from the clustered BAM

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 04_CLUSTERING
job_header "cluster stats"
activate_env
need_tools samtools
mkdir -p results/qc

for G in ALL $(sheet_groups); do
    REPORT="results/${G}/clustered_${G}.cluster_report.csv"
    BAM="results/${G}/clustered_${G}.bam"
    OUT="results/qc/cluster_stats_${G}.tsv"
    if [ ! -f "$REPORT" ] || [ ! -f "$BAM" ]; then
        echo "$G: not clustered, skipping"
        continue
    fi
    echo "$G"
    : > "$OUT"

    # FLNC input = sum of poly(A) FLNC counts from each refine report
    FLNC=0
    while read -r f; do
        j="${f%.bam}.filter_summary.report.json"
        v=$(awk -F': ' '/"num_reads_flnc_polya"/{getline; getline; gsub(/[^0-9]/,"",$2); print $2}' "$j" 2>/dev/null)
        FLNC=$((FLNC + ${v:-0}))
    done < "${FOFN_DIR}/flnc_${G}.fofn"
    printf "flnc_input\t%s\n" "$FLNC" >> "$OUT"

    # reads per cluster (skip header)
    awk -F, 'NR > 1 { c[$1]++ }
    END {
        for (k in c) { n = c[k]; total++
            if (n == 1) s1++
            else if (n <= 5) s2++
            else if (n <= 10) s3++
            else if (n <= 50) s4++
            else s5++
        }
        printf "total_isoforms\t%d\n", total
        printf "support_1\t%d\n", s1
        printf "support_2_5\t%d\n", s2
        printf "support_6_10\t%d\n", s3
        printf "support_11_50\t%d\n", s4
        printf "support_50plus\t%d\n", s5
    }' "$REPORT" >> "$OUT"

    # transcript length
    samtools view "$BAM" | awk '{
        L = length($10); n++; sum += L
        if (n == 1 || L < mn) mn = L
        if (L > mx) mx = L
        if (L <= 500) b1++; else if (L <= 1000) b2++; else if (L <= 2000) b3++
        else if (L <= 3000) b4++; else if (L <= 5000) b5++; else if (L <= 10000) b6++; else b7++
    }
    END {
        printf "len_n\t%d\n", n
        printf "mean_length\t%.2f\n", (n ? sum / n : 0)
        printf "min_length\t%d\n", mn
        printf "max_length\t%d\n", mx
        printf "len_0_500\t%d\n", b1
        printf "len_501_1000\t%d\n", b2
        printf "len_1001_2000\t%d\n", b3
        printf "len_2001_3000\t%d\n", b4
        printf "len_3001_5000\t%d\n", b5
        printf "len_5001_10000\t%d\n", b6
        printf "len_10000plus\t%d\n", b7
    }' >> "$OUT"

    sed 's/^/    /' "$OUT"
done
echo "next: python CLUSTER_QC.py"
job_footer
