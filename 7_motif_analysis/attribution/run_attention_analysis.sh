#!/bin/bash
#SBATCH --job-name=attn_attr
#SBATCH --ntasks=1
#SBATCH --output=attn_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=01:00:00

# ProstT5 self-attention vs classifier IG, per residue (showcase proteins).
# Needs ig_output_logit/ (run run_integrated_gradients.sh first).
cd /home/nroethler/Code/motif_analysis/attribution
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python -u attention_analysis.py --fasta showcase.fasta --out-dir attn_output
echo "attention analysis done -- overlays + attn_vs_ig_summary.csv in attn_output/"
