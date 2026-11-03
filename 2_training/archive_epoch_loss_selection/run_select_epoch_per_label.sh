#!/bin/bash
#SBATCH --job-name=select_epoch_label
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=select_epoch_per_label.log

# Choose the epoch of every (embedding, label) from the BCE validation loss alone (CPU, about a minute); no MCC, no loss selection.
#   sbatch run_select_epoch_per_label.sh
#   SELECT_ARGS="--dir <sweep>/train_val_loss_output --out <dir> --splits 6" sbatch run_select_epoch_per_label.sh     # select on the first 6 splits only; the bake-off keeps 7-16 unseen
# Reads train_val_loss_output/ (train_val_loss_curves.py), writes epoch_selection_per_label_output/ (selection_per_label_<set>.json is the bake-off's LOSS_MAP_JSON).
# Environment: the conda environment that is already active, else the one `make env` builds (SMSTPRED_ENV renames it).
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
python -u test_select_epoch_per_label.py
python -u select_epoch_per_label.py --set "${SELECT_SET:-sam_full}" ${SELECT_ARGS:-}
