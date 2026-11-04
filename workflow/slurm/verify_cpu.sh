#!/bin/bash
#SBATCH --job-name=verify_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=verify_cpu_%j.log

# `make verify` as one job: runs every CPU stage on the release assembled from the archived record (make data) and compares with the archived results. Submit from the repository root:
#   sbatch workflow/slurm/verify_cpu.sh        (SMSTPRED_ENV=<conda env> to use another environment than `smstpred`)
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
[ -f Makefile ] || { echo "submit from the repository root" >&2; exit 1; }
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
make verify
