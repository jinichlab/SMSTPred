#!/bin/bash
#SBATCH --job-name=render_dna_motifs
#SBATCH --ntasks=1
#SBATCH --output=render_dna_motifs.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00

# Render one ChimeraX figure per DNA STREME motif (offscreen, software render,
# no GPU). Motif sidechains as sticks with N=blue O=red S=yellow heteroatoms.
cd /home/nroethler/Code/motif_analysis/structures/DNA
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python render_dna_motifs.py
echo "DNA motif figures done -- structures/DNA/motif*/"
