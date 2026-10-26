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
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python -u train_val_loss_curves.py --task-id $SLURM_ARRAY_TASK_ID
