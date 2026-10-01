# lib/functions.sh
# Helpers shared by all step scripts. Sources config.sh for you.
# You should not need to edit this file; change config.sh instead.
#
# Each step script starts with:
#   HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
#   source "${HERE}/../lib/functions.sh" || exit 1
# SLURM copies the script to a spool dir, so $0 is useless inside a job.
# That is why you have to sbatch from inside the step folder.

_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_LIB_DIR}/../config.sh" || { echo "ERROR: could not source config.sh"; exit 1; }

die()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }

# relative paths in config.sh are relative to PROJECT_DIR
abspath() {
    local p="$1"
    [ -z "$p" ] && { echo ""; return; }
    case "$p" in /*) echo "$p" ;; *) echo "${PROJECT_DIR}/$p" ;; esac
}

SAMPLE_SHEET="$(abspath "$SAMPLE_SHEET")"
SKERA_ADAPTERS="$(abspath "$SKERA_ADAPTERS")"
LIMA_PRIMERS="$(abspath "$LIMA_PRIMERS")"
GENOME_FA="$(abspath "$GENOME_FA")"
ANNOTATION_GTF="$(abspath "$ANNOTATION_GTF")"
CAGE_BED="$(abspath "$CAGE_BED")"
POLYA_MOTIFS="$(abspath "$POLYA_MOTIFS")"
for _i in "${!HIFI_BAMS[@]}"; do HIFI_BAMS[$_i]="$(abspath "${HIFI_BAMS[$_i]}")"; done
unset _i

# pbmm2 index lives next to the genome
GENOME_MMI="${GENOME_FA%.*}.ISOSEQ.mmi"
# pigeon prepare turns foo.gtf into foo.sorted.gtf (+ .pgi)
SORTED_GTF="${ANNOTATION_GTF%.gtf}.sorted.gtf"
# same idea for the optional CAGE bed
CAGE_SORTED=""
[ -n "$CAGE_BED" ] && CAGE_SORTED="${CAGE_BED%.bed}.sorted.bed"

# step folders
SEG_DIR="${PROJECT_DIR}/01_SEGMENTATION"
LIMA_DIR="${PROJECT_DIR}/02_LIMA_REFINE"
FOFN_DIR="${PROJECT_DIR}/03_MERGE_FOFN/fofn_files"
CLUST_DIR="${PROJECT_DIR}/04_CLUSTERING"
MAP_DIR="${PROJECT_DIR}/05_MAPPING"
COLL_DIR="${PROJECT_DIR}/06_COLLAPSE"
CLASS_DIR="${PROJECT_DIR}/07_CLASSIFY"

# make sure the job was submitted from the right folder, otherwise the
# relative --output=logs/... paths and results/ end up somewhere odd
in_step_dir() {
    local want="$1"
    [ "$(basename "$HERE")" = "$want" ] || \
        die "submit this script from inside ${want}/ (you are in ${HERE})"
    cd "$HERE" || die "cannot cd to $HERE"
}

activate_env() {
    # conda activate trips over set -u, so relax it for a moment
    local had_u=0
    case "$-" in *u*) had_u=1; set +u ;; esac
    if [ -n "${MODULES}" ]; then
        module load ${MODULES} || warn "module load ${MODULES} failed"
    fi
    [ -f "$CONDA_SH" ] || die "CONDA_SH not found: $CONDA_SH (edit config.sh)"
    source "$CONDA_SH"
    conda activate "$CONDA_ENV" || die "could not activate conda env '$CONDA_ENV'"
    [ $had_u -eq 1 ] && set -u
    return 0
}

need_tools() {
    local t
    for t in "$@"; do
        command -v "$t" >/dev/null 2>&1 || die "'$t' not on PATH. Is the env right? (CONDA_ENV=$CONDA_ENV)"
    done
}

need_file() {
    # need_file <path> [hint]
    [ -s "$1" ] || die "missing or empty: $1${2:+  ($2)}"
}

job_header() {
    echo "=================================================="
    echo "$1"
    echo "job ${SLURM_JOB_ID:-local}  task ${SLURM_ARRAY_TASK_ID:-na}  node ${SLURMD_NODENAME:-$(hostname)}"
    echo "started $(date)"
    echo "=================================================="
}

job_footer() { echo "finished $(date)"; }

# pick the Nth (1-based) entry of a list for array jobs
# usage: X=$(array_pick "${LIST[@]}")
array_pick() {
    local n="${SLURM_ARRAY_TASK_ID:-}"
    [ -n "$n" ] || die "this script is meant to run as an array job (sbatch --array=1-N)"
    local list=("$@")
    [ "$n" -ge 1 ] && [ "$n" -le "${#list[@]}" ] || \
        die "array task $n is out of range (only ${#list[@]} inputs). Fix --array."
    echo "${list[$((n-1))]}"
}

# m84xxx_s1.hifi_reads.bcM0001.bam -> m84xxx_s1.bcM0001
cell_id() {
    local b
    b="$(basename "$1" .bam)"
    b="${b/.hifi_reads/}"
    b="${b%.segmented}"
    echo "$b"
}

# ---------------------------------------------------------------------------
# sample sheet
# ---------------------------------------------------------------------------
# rows only, no comments/header, CR stripped. cols: barcode sample group cell
sheet_rows() {
    [ -s "$SAMPLE_SHEET" ] || die "sample sheet not found: $SAMPLE_SHEET"
    awk -F'\t' 'BEGIN{OFS="\t"} {sub(/\r$/,"")} /^#/ || NF<3 || $1=="barcode" {next}
                { c=(NF>=4 && $4!="") ? $4 : "all"; print $1,$2,$3,c }' "$SAMPLE_SHEET"
}
sheet_samples() { sheet_rows | awk -F'\t' '!s[$2]++{print $2}'; }
sheet_groups()  { sheet_rows | awk -F'\t' '!s[$3]++{print $3}'; }
group_of()      { sheet_rows | awk -F'\t' -v s="$1" '$2==s{print $3; exit}'; }

# sample name for a given cell id + barcode pair (empty if not in sheet)
sample_for() {
    sheet_rows | awk -F'\t' -v cell="$1" -v pair="$2" \
        '$1==pair && ($4=="all" || index(cell,$4)>0) {print $2; exit}'
}

check_sheet() {
    local n bad
    n=$(sheet_rows | wc -l)
    [ "$n" -gt 0 ] || die "no samples in $SAMPLE_SHEET"
    bad=$(sheet_rows | awk -F'\t' '$2 ~ /[ ,]/ {print $2}')
    [ -z "$bad" ] || die "sample names cannot have spaces or commas: $bad"
    bad=$(sheet_rows | awk -F'\t' '$3=="ALL"{print $2}')
    [ -z "$bad" ] || die "'ALL' is reserved, pick another group name (samples: $bad)"
    bad=$(sheet_rows | awk -F'\t' 'seen[$1 FS $4]++ {print $1" ("$4")"}')
    [ -z "$bad" ] || die "barcode listed twice for the same cell: $bad"
}

# GROUP argument is ALL or one of the groups in the sheet
check_group() {
    local g="$1"
    [ -n "$g" ] || die "no group given (use ALL or one of: $(sheet_groups | tr '\n' ' '))"
    [ "$g" = "ALL" ] && return 0
    sheet_groups | grep -qx -- "$g" || die "group '$g' is not in the sample sheet"
}

# lima writes <cell>.<pair>.bam ; the pair has no dots, the cell id can
demux_bams() {
    ls "${LIMA_DIR}"/results/demux/*.bam 2>/dev/null | grep -v -e unbarcoded -e '\.flnc\.' | sort
}
split_demux_name() {
    # prints "<cell>\t<pair>" for .../<cell>.<pair>[.flnc].bam
    local b
    b="$(basename "$1")"; b="${b%.bam}"; b="${b%.flnc}"
    printf "%s\t%s\n" "${b%.*}" "${b##*.}"
}

count_reads() { samtools view -c -@ 4 "$1" 2>/dev/null; }
