#!/bin/bash
#SBATCH --job-name=collapse_stats
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 06b (QC): summarise the collapse GFF + flnc_count.txt for COLLAPSE_QC.py
#
#   cd 06_COLLAPSE
#   sbatch 06b_EXTRACT_COLLAPSE_STATS.sh      (after 05d has run)
#   (when it finishes)  python COLLAPSE_QC.py
#
# Writes results/qc/<BUILD_LABEL>/collapse_stats_<GROUP>.tsv
#   pre_isoforms   mapped transcripts before collapse (from 05d)
#   post_isoforms, loci, ipl_* (isoforms per locus), span_*, exon_mono/multi
#   supp_*         FLNC reads per isoform
#   ps_<sample>_flnc / ps_<sample>_iso   per-sample reads / isoforms seen

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 06_COLLAPSE
job_header "collapse stats: $BUILD_LABEL"

QC="results/qc/${BUILD_LABEL}"
MAPQC="${MAP_DIR}/results/qc/${BUILD_LABEL}"
mkdir -p "$QC"

for G in ALL $(sheet_groups); do
    base="results/${BUILD_LABEL}/${G}/clustered_${G}.${BUILD_LABEL}.collapsed"
    GFF="${base}.gff"; FLNC="${base}.flnc_count.txt"
    OUT="${QC}/collapse_stats_${G}.tsv"
    [ -f "$GFF" ] || { echo "$G: not collapsed, skipping"; continue; }
    : > "$OUT"

    PRE=$(awk -F'\t' '$1=="flagstat_primary_mapped"{print $2; exit}' "${MAPQC}/mapstats_${G}.tsv" 2>/dev/null)
    [ -n "$PRE" ] || warn "no mapstats for $G, run 05d first (pre_isoforms set to 0)"
    printf "pre_isoforms\t%s\n" "${PRE:-0}" >> "$OUT"

    # GFF lists each transcript line followed by its exons
    awk -F'\t' '
      $3 == "transcript" {
          if (have) { if (ec == 1) mono++; else if (ec > 1) multi++ }
          niso++
          tid = $9; sub(/.*transcript_id "/, "", tid); sub(/".*/, "", tid)   # PB.X.Y
          n = split(tid, p, "."); if (n >= 2) loc[p[2]]++
          span = $5 - $4 + 1; sumspan += span
          if (span <= 200) s0++; else if (span <= 500) s1++; else if (span <= 1000) s2++
          else if (span <= 2000) s3++; else if (span <= 5000) s4++; else if (span <= 10000) s5++
          else if (span <= 15000) s6++; else s7++
          ec = 0; have = 1
      }
      $3 == "exon" { ec++ }
      END {
          if (have) { if (ec == 1) mono++; else if (ec > 1) multi++ }
          for (k in loc) { nloc++; c = loc[k]
              if (c == 1) i1++; else if (c == 2) i2++; else if (c == 3) i3++
              else if (c <= 5) i45++; else if (c <= 10) i610++; else if (c <= 20) i1120++; else i20++ }
          printf "post_isoforms\t%d\nloci\t%d\n", niso, nloc
          printf "ipl_1\t%d\nipl_2\t%d\nipl_3\t%d\nipl_4_5\t%d\n", i1, i2, i3, i45
          printf "ipl_6_10\t%d\nipl_11_20\t%d\nipl_20plus\t%d\n", i610, i1120, i20
          printf "span_le200\t%d\nspan_201_500\t%d\nspan_501_1000\t%d\n", s0, s1, s2
          printf "span_1001_2000\t%d\nspan_2001_5000\t%d\nspan_5001_10k\t%d\n", s3, s4, s5
          printf "span_10k_15k\t%d\nspan_15kplus\t%d\n", s6, s7
          printf "span_mean\t%.1f\n", (niso ? sumspan / niso : 0)
          printf "exon_mono\t%d\nexon_multi\t%d\n", mono, multi
      }' "$GFF" >> "$OUT"

    if [ -f "$FLNC" ]; then
        # comma separated in isoseq 4.x; split on tab too just in case
        awk -F'[,\t]' '
          /^#/ { next }
          !hdr { ncol = NF; for (c = 2; c <= NF; c++) name[c] = $c; hdr = 1; next }
          {
              tot = 0
              for (c = 2; c <= NF; c++) { v = $c + 0; tot += v; colsum[c] += v; if (v > 0) colnz[c]++ }
              if (tot <= 1) t1++; else if (tot <= 5) t25++; else if (tot <= 10) t610++
              else if (tot <= 50) t1150++; else if (tot <= 100) t51100++; else t100++
          }
          END {
              printf "supp_1\t%d\nsupp_2_5\t%d\nsupp_6_10\t%d\n", t1, t25, t610
              printf "supp_11_50\t%d\nsupp_51_100\t%d\nsupp_100plus\t%d\n", t1150, t51100, t100
              for (c = 2; c <= ncol; c++) {
                  printf "ps_%s_flnc\t%d\n", name[c], colsum[c]
                  printf "ps_%s_iso\t%d\n", name[c], colnz[c]
              }
          }' "$FLNC" >> "$OUT"
    else
        warn "$FLNC missing, no support/per-sample numbers for $G"
    fi
    echo "$G: $(grep -E '^(pre|post)_isoforms' "$OUT" | tr '\n' ' ')"
done
echo "next: python COLLAPSE_QC.py"
job_footer
