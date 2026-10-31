#!/bin/bash
#SBATCH --job-name=verify_gpu
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=verify_gpu_%j.log

# `make verify-gpu` as one job: re-trains the deployment model and the MLP bake-off, re-predicts all proteomes, re-runs the loss-curve sweeps
# (as a SLURM array, ten sweeps in parallel) and compares the outcomes with the stated tolerances. Needs the data in place (`make data`) and the
# repository environment active (or `make env` done). Submit from the repository root:   sbatch workflow/slurm/verify_gpu.sh
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
[ -f Makefile ] || { echo "submit from the repository root" >&2; exit 1; }
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
make verify-gpu
