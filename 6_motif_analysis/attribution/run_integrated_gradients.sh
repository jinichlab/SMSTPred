#!/bin/bash
#SBATCH --job-name=ig_attr
#SBATCH --ntasks=1
#SBATCH --output=ig_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=02:00:00

# Per-residue Integrated-Gradients attribution for the ProstT5 -> MLP classifier.
#   1) build showcase FASTA from the labelled CSV
#   2) ProstT5 per-residue embeddings (prott5 env, GPU)  <- /residue groups needed
#   3) IG attribution + per-residue tracks (classifier_env)
set -e
cd /home/nroethler/Code/motif_analysis/attribution
source ~/miniconda3/etc/profile.d/conda.sh

# 1) FASTA of the curated showcase (one+ per substrate class, known catalytic motif)
conda activate classifier_env
python prep_showcase_fasta.py --out showcase.fasta

# 2) per-residue ProstT5 embeddings (per-residue is ON by default in embed_prostt5.py)
conda activate prott5
python /home/nroethler/Code/prostt5/embed_prostt5.py \
    --fasta showcase.fasta \
    --out showcase_prostt5.h5

# 3) Integrated Gradients (mean-MTase baseline = "what makes this a <class>-MTase").
#    logit space: unsaturated gradients that work for every class incl. the
#    majority (RNA) -- prob space saturates at high confidence + the mean baseline
#    is degenerate for RNA. Pass --space prob to attribute the deployed score.
conda activate classifier_env
python integrated_gradients.py \
    --residue-h5 showcase_prostt5.h5 \
    --fasta showcase.fasta \
    --label auto \
    --baseline mean \
    --space logit \
    --steps 64 \
    --out-dir ig_output

echo "IG done -- tracks + ig_summary.csv in ig_output/"
