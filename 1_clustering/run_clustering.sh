#!/bin/bash
#SBATCH --job-name=clustering
#SBATCH --ntasks=1
#SBATCH --output=run_clustering_output.log
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python sam_clustering.py
python up_clustering.py
