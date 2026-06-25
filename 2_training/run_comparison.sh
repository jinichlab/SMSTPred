#!/bin/bash
#SBATCH --job-name=sam_clf_cmp
#SBATCH --ntasks=1
#SBATCH --output=run_comparison_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python plot_cluster_mcc_zoom.py
