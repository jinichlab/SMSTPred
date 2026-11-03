#!/bin/bash
#SBATCH --job-name=select_epoch_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=select_epoch_loss.log

# Choose the training epoch and the loss of the MLP classifiers from the sweep that run_train_val_loss.sh recorded (CPU, about a minute).
#   sbatch run_select_epoch_loss.sh                                   # the SAM sweep (sam_full): epoch from the validation-loss minima, loss per embedding and label by MCC there
#   SELECT_SET=up sbatch run_select_epoch_loss.sh                     # the UP ProstT5 sweep
#   SELECT_ARGS="--loss-scope hierarchical" sbatch run_select_epoch_loss.sh   # one global loss, changed only where a cell's gain is real (see select_epoch_and_loss.py --help)
#   SELECT_ARGS="--epoch 300" ...                                     # fix the epoch by hand
# Reads train_val_loss_output/ (written by train_val_loss_curves.py: <set>__<embedding>__curves.csv.gz and __ckpt_mcc.csv.gz), writes epoch_loss_selection_output/.
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
python -u test_select_epoch_and_loss.py
python -u select_epoch_and_loss.py --set "${SELECT_SET:-sam_full}" ${SELECT_ARGS:-}
