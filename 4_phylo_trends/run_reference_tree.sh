#!/bin/bash
#SBATCH --job-name=ref_tree
#SBATCH --ntasks=1
#SBATCH --output=reference_tree_output.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python reference_tree.py
