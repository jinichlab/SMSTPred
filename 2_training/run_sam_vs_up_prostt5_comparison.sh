#!/bin/bash
#SBATCH --job-name=sam_vs_up
#SBATCH --ntasks=1
#SBATCH --output=sam_vs_up_prostt5_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=04:00:00

# SAM-trained vs UP-trained ProstT5 MLP on a held-out SAM test set (5-fold
# cluster-disjoint). MVP: only the exact test entries are removed from UP training.
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python -u sam_vs_up_prostt5_comparison.py
