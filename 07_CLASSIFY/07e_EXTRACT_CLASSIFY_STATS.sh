#!/bin/bash
#SBATCH --job-name=classify_stats
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 07e (QC): parse the pigeon summaries for CLASSIFY_QC.py
#
#   cd 07_CLASSIFY
#   sbatch 07e_EXTRACT_CLASSIFY_STATS.sh      (or just bash it, it is quick)
#   python CLASSIFY_QC.py
#
# Writes results/qc/<BUILD_LABEL>/classify_stats_<GROUP>.tsv
#   pre_*    from <prefix>.summary.txt (before filter)
#   post_*   from *_classification.filtered.summary.txt
#   reason_* filter reasons
#   ps_<sample>_flpass / ps_<sample>_isopass   per-sample reads / isoforms that passed
# Category keys: FSM ISM NIC NNC Antisense GenicIntron GenicGenomic Intergenic Other
# (<pre|post>_<CAT> by isoform, <pre|post>read_<CAT> by read)

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 07_CLASSIFY
job_header "classify stats: $BUILD_LABEL"

QC="results/qc/${BUILD_LABEL}"
mkdir -p "$QC"

# "Label : 123 (4.5%)" -> 123
num() {
    awk -F: -v L="$2" '{ k = $1; gsub(/^[ \t]+|[ \t]+$/, "", k)
        if (k == L) { v = $2; sub(/\(.*/, "", v); gsub(/[, ]/, "", v); print v + 0; exit } }' "$1"
}

# categories, by isoform and by read
cats() {
    awk -F: -v P="$2" '
    function key(s) { gsub(/^[ \t]+|[ \t]+$/, "", s)
        if (s == "Full splice match") return "FSM"
        if (s == "Incomplete splice match") return "ISM"
        if (s == "Novel in catalog") return "NIC"
        if (s == "Novel not in catalog") return "NNC"
        if (s == "Antisense") return "Antisense"
        if (s == "Genic intron") return "GenicIntron"
        if (s == "Genic genomic") return "GenicGenomic"
        if (s == "Intergenic") return "Intergenic"
        if (s == "Other") return "Other"
        return "" }
    /Classifications, by isoform/ { sec = P; next }
    /Classifications, by read/    { sec = P "read"; next }
    /^(Junctions|Genes|Filter reasons|RT Switching)/ { sec = "" }
    sec != "" && /:/ { k = key($1); if (k == "") next
        v = $2; sub(/\(.*/, "", v); gsub(/[, ]/, "", v)
        printf "%s_%s\t%d\n", sec, k, v }' "$1"
}

# Known/Novel genes under the "Genes" heading
genes() {
    awk -v P="$2" '/^Genes/ { g = 1; next }
        g && /^Known/ { v = $0; sub(/.*:[ ]*/, "", v); sub(/\(.*/, "", v); gsub(/[, ]/, "", v); printf "%s_genes_known\t%d\n", P, v }
        g && /^Novel/ { v = $0; sub(/.*:[ ]*/, "", v); sub(/\(.*/, "", v); gsub(/[, ]/, "", v); printf "%s_genes_novel\t%d\n", P, v; exit }' "$1"
}

for G in ALL $(sheet_groups); do
    D="results/${BUILD_LABEL}/${G}"
    CS="${D}/clustered_${G}.summary.txt"
    FS="${D}/clustered_${G}_classification.filtered.summary.txt"
    FC="${D}/clustered_${G}_classification.filtered_lite_classification.txt"
    OUT="${QC}/classify_stats_${G}.tsv"
    [ -f "$CS" ] || { echo "$G: not classified, skipping"; continue; }
    : > "$OUT"

    printf "pre_input\t%s\n" "$(num "$CS" 'Input transcripts')" >> "$OUT"
    cats "$CS" pre >> "$OUT"
    genes "$CS" pre >> "$OUT"

    if [ -f "$FS" ]; then
        {
            printf "post_input\t%s\n"   "$(num "$FS" 'Input')"
            printf "post_passed\t%s\n"  "$(num "$FS" 'Passed')"
            printf "post_removed\t%s\n" "$(num "$FS" 'Removed')"
            printf "reason_intrapriming\t%s\n" "$(num "$FS" 'Intrapriming')"
            printf "reason_monoexonic\t%s\n"   "$(num "$FS" 'Monoexonic')"
            printf "reason_rtswitch\t%s\n"     "$(num "$FS" 'RT switching')"
            printf "reason_lowcov\t%s\n"       "$(num "$FS" 'Low coverage/non-canonical')"
        } >> "$OUT"
        cats "$FS" post >> "$OUT"
        genes "$FS" post >> "$OUT"
    else
        warn "no filter summary for $G"
    fi

    # per-sample FL.<sample> columns in the filtered classification
    if [ -f "$FC" ]; then
        awk -F'\t' '
            NR == 1 { for (i = 1; i <= NF; i++) if ($i ~ /^FL\./) { idx[++n] = i; name[i] = substr($i, 4) }; next }
            { for (j = 1; j <= n; j++) { c = idx[j]; v = $c + 0; fl[c] += v; if (v > 0) iso[c]++ } }
            END { for (j = 1; j <= n; j++) { c = idx[j]
                    printf "ps_%s_flpass\t%d\n", name[c], fl[c]
                    printf "ps_%s_isopass\t%d\n", name[c], iso[c] } }' "$FC" >> "$OUT"
    fi
    echo "$G: input $(num "$CS" 'Input transcripts'), passed $( [ -f "$FS" ] && num "$FS" 'Passed')"
done
echo "next: python CLASSIFY_QC.py"
job_footer
