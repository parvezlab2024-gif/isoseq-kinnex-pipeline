#!/bin/bash
#SBATCH --job-name=pigeon_report
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 07d: pigeon report -> saturation table (genes / isoforms found as a
# function of reads). Tells you whether more sequencing would find much more.
#
#   cd 07_CLASSIFY
#   sbatch 07d_PIGEON_REPORT.sh ALL
#
# Uses --exclude-singletons, since singletons were kept at clustering and
# would otherwise make the curve look like it never flattens.
# CLASSIFY_QC.py plots the table if it is there.

set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${HERE}/../lib/functions.sh" || exit 1
in_step_dir 07_CLASSIFY

GROUP="${1:-ALL}"
check_group "$GROUP"
job_header "pigeon report: $GROUP ($BUILD_LABEL)"
activate_env
need_tools pigeon

OUTDIR="results/${BUILD_LABEL}/${GROUP}"
FILT="${OUTDIR}/clustered_${GROUP}_classification.filtered_lite_classification.txt"
OUT="${OUTDIR}/clustered_${GROUP}.saturation.txt"
need_file "$FILT" "run 07c $GROUP first"

pigeon report --exclude-singletons "$FILT" "$OUT" || die "pigeon report failed"
need_file "$OUT"
head -3 "$OUT"; echo "..."; tail -2 "$OUT"
job_footer
