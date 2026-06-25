#!/bin/bash
#SBATCH --job-name=up_prostt5_loss
#SBATCH --ntasks=1
#SBATCH --output=up_prostt5_loss_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# UP (all-UniProt) ProstT5 MLP loss/epoch comparison, cluster-disjoint (6 iters).
# Requires files/prostt5_UP_swiss.h5 (run prostt5/run_prostt5_up.sh first).
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python -u mlp_loss_comparison_prostt5.py
