#!/bin/bash
#SBATCH --job-name=tv_loss
#SBATCH --ntasks=1
#SBATCH --array=0-5,12-15%6
#SBATCH --output=train_val_loss_%A_%a.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00

# Train vs. validation loss per epoch for the S2/S5 sweeps (overfitting check; Figs S3 and S6).
# Array: 0-5 SAM full x 6 embeddings, 12-15 UP ProstT5 x 4 losses (`python train_val_loss_curves.py --list`).
# Tasks 6-11 (SAM balanced subset, not used in the manuscript) are left out; add them to --array to run them.
# Inputs are read from <repo>/data, or from $SMSTPRED_DATA.
# Resumable: finished (iteration, loss) pairs are skipped on resubmit.
# Afterwards (CPU): python plot_train_val_loss_curves.py
# Environment: the conda environment that is already active, else the one `make env` builds (SMSTPRED_ENV renames it).
# Runs as `sbatch run_*.sh` (the job inherits the submitting shell's environment) or as plain `bash run_*.sh`.
# The scripts use relative paths, so start in this script's folder (under sbatch $0 is a spooled copy: ask SLURM where it came from).
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi

# One array task per sweep under sbatch; without SLURM_ARRAY_TASK_ID (plain bash) the ten sweeps run one after the other.
for task in ${SLURM_ARRAY_TASK_ID:-0 1 2 3 4 5 12 13 14 15}; do
    python -u train_val_loss_curves.py --task-id "$task"
done
