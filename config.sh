# config.sh
# Settings for the Iso-Seq pipeline. Every step script sources this file,
# so this (plus samples.tsv) is normally the only thing you need to edit.
# Written by Talha Lone (Parvez lab). Last updated 2026-09-17.
#
# Paths can be absolute or relative to PROJECT_DIR.
# Nothing in here runs anything, it only sets variables.

# ---------------------------------------------------------------------------
# 1. Where the pipeline lives
# ---------------------------------------------------------------------------
# Leave this as-is unless you move config.sh somewhere else. It resolves to
# the folder that holds config.sh.
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# 2. Software
# ---------------------------------------------------------------------------
# conda.sh from your mamba/conda install (this is the Quest path).
CONDA_SH="/hpc/software/mamba/23.1.0/etc/profile.d/conda.sh"
# Name (or full path) of the env built from SETUP/isoseq_environment.yml
CONDA_ENV="pacbio_workflow"
# Modules to load before each job. Leave empty if samtools is in the env
# (it is, if you built it from our yml).
MODULES=""

# ---------------------------------------------------------------------------
# 3. Raw data
# ---------------------------------------------------------------------------
# One HiFi BAM per SMRT cell (or per Kinnex/MAS-Seq index if the core split
# them, e.g. *.hifi_reads.bcM0001.bam). Add one line per file.
# 01a/01b/02a run as SLURM arrays, one task per entry here.
HIFI_BAMS=(
    "RAW/CHANGE_ME/hifi_reads/m84XXX_XXXXXX_XXXXXX_s1.hifi_reads.bcM0001.bam"
)

# Sample sheet: barcode pair -> sample name -> group. See samples.tsv.
SAMPLE_SHEET="samples.tsv"

# ---------------------------------------------------------------------------
# 4. Library prep
# ---------------------------------------------------------------------------
# Skera adapters. mas8 = Kinnex full-length RNA (8-fold). If your library is
# MAS-Seq 16-fold, download mas16_primers.fasta from PacBio and point here.
SKERA_ADAPTERS="01_SEGMENTATION/adapters/mas8_primers.fasta"
EXPECTED_ARRAY_SIZE=8          # used by SKERA_QC.py only

# cDNA barcoded primers (Iso-Seq v2 / Kinnex, bc01-bc12 + IsoSeqX_3p).
LIMA_PRIMERS="02_LIMA_REFINE/primers/IsoSeq_v2_primers_12.fasta"

# isoseq refine
MIN_POLYA_LENGTH=20

# ---------------------------------------------------------------------------
# 5. Reference (see REFERENCE/README.md before filling these in)
# ---------------------------------------------------------------------------
# Short name for this genome + annotation. Used in output folder/file names,
# so no spaces. e.g. GRCh38_v47, GRCm39_vM38, GRCz11_r113
BUILD_LABEL="CHANGE_ME"
GENOME_FA="REFERENCE/CHANGE_ME.fa"
# Pigeon-compatible GTF (NOT the .sorted.gtf; 07b makes that for you).
ANNOTATION_GTF="REFERENCE/CHANGE_ME.gtf"

# Optional extras for pigeon classify. Leave empty to skip.
CAGE_BED=""          # CAGE / refTSS peaks, same assembly + chrom names as GENOME_FA
POLYA_MOTIFS=""      # ranked polyA motif list, one per line

# Regex for the "main" chromosomes in the mapping QC plot. Anything that does
# not match is lumped into "other". Default handles chr1/1, X, Y, M/MT.
# NCBI RefSeq genomes use accessions (NC_000001.11), which won't match;
# for those use something like '^NC_'.
MAIN_CHROM_REGEX='^(chr)?([0-9]+|[XYZW]|M|MT)$'

# ---------------------------------------------------------------------------
# 6. Tool options
# ---------------------------------------------------------------------------
# --singletons keeps single-read clusters. They are dropped later by
# pigeon report --exclude-singletons and by any count filter downstream.
CLUSTER_OPTS="--singletons"
# We keep 5' truncated isoforms separate (did not merge them into the longer
# model). Remove the flag to use the isoseq default.
COLLAPSE_OPTS="--do-not-collapse-extra-5exons"

# Threads. Keep these in line with --cpus-per-task in each script.
THREADS_SKERA=8
THREADS_LIMA=8
THREADS_REFINE=4
THREADS_CLUSTER=16
THREADS_MAP=16
THREADS_COLLAPSE=8
THREADS_PIGEON=8
