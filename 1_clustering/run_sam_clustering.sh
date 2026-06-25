#!/bin/bash
#SBATCH --job-name=sam_reclust
#SBATCH --ntasks=1
#SBATCH --output=sam_reclust_output.log
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python sam_clustering.py
