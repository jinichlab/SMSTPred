#!/bin/bash
#SBATCH --job-name=domain_enrich
#SBATCH --ntasks=1
#SBATCH --output=domain_enrichment_output.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python 3_domain_enrichment.py
