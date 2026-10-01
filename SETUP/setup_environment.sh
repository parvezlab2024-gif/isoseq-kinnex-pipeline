#!/bin/bash
# Build the conda env from isoseq_environment.yml.
#   bash SETUP/setup_environment.sh              -> env named pacbio_workflow
#   bash SETUP/setup_environment.sh my_env_name  -> different name
# If you change the name, change CONDA_ENV in config.sh to match.
# Takes a while; run it in an interactive job rather than on a login node
# if your cluster is strict about that.

set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ENV_FILE="${HERE}/isoseq_environment.yml"
ENV_NAME="${1:-pacbio_workflow}"

if command -v mamba >/dev/null 2>&1; then
    CONDA_CMD=mamba
elif command -v conda >/dev/null 2>&1; then
    CONDA_CMD=conda
else
    echo "no mamba or conda on PATH. On Quest: module load mamba (or source its conda.sh)"
    exit 1
fi
echo "using $CONDA_CMD, env name $ENV_NAME"

if $CONDA_CMD env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "env $ENV_NAME already exists."
    echo "to update it:   $CONDA_CMD env update -n $ENV_NAME -f $ENV_FILE"
    echo "to start over:  $CONDA_CMD env remove -n $ENV_NAME   and rerun this"
    exit 0
fi

$CONDA_CMD env create -n "$ENV_NAME" -f "$ENV_FILE"

# quick check that everything landed
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"
echo ""
for t in skera lima isoseq pbmm2 pigeon; do
    printf "%-8s %s\n" "$t" "$($t --version 2>/dev/null | head -1)"
done
printf "%-8s %s\n" samtools "$(samtools --version | head -1)"
python -c "import matplotlib, pandas, numpy; print('python plotting libs ok')"
echo ""
echo "done. activate with: conda activate $ENV_NAME"
