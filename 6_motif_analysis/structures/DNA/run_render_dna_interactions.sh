#!/bin/bash
#SBATCH --job-name=render_dna_interact
#SBATCH --ntasks=1
#SBATCH --output=render_dna_interactions.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00

# Interaction-focused per-motif figures: motif sticks + only its binding partner
# (DNA / cofactor) + H-bond dashes, bulk protein hidden. ChimeraX offscreen, no GPU.
cd /home/nroethler/Code/motif_analysis/structures/DNA
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python render_dna_motif_interactions.py
echo "DNA motif interaction figures done -- structures/DNA/motif*/*_interaction.png"
