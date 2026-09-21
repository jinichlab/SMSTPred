#!/bin/bash
#SBATCH --job-name=verify_upgrade
#SBATCH --ntasks=1
#SBATCH --output=verify_upgrade_motifs_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# Confirm each apo/AF-motif upgrade structure's protein actually carries the target
# motif: fetch its UniProt sequence and FIMO-scan with the class's STREME motifs.
# Needs outbound HTTPS (rest.uniprot.org) + fimo (MEME suite, classifier_env).
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python verify_upgrade_motifs.py
