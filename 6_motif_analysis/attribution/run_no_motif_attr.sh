#!/bin/bash
#SBATCH --job-name=nomotif_ig
#SBATCH --ntasks=1
#SBATCH --output=no_motif_attr_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=02:00:00

# Per-residue IG over the 3,023 no-motif reps + trend aggregation.
#   1) FASTA + manifest (no_motif_prep.py)
#   2) ProstT5 per-residue embeddings for just these reps (~1 min)
#   3) IG toward each rep's class + AA-enrichment / peakedness / top-window trends
cd /home/nroethler/Code/motif_analysis/attribution
source ~/miniconda3/etc/profile.d/conda.sh

conda activate classifier_env
python no_motif_prep.py

conda activate prott5
python /home/nroethler/Code/prostt5/embed_prostt5.py \
    --fasta no_motif.fasta --out no_motif_prostt5.h5     # per-residue ON by default

conda activate classifier_env
python -u no_motif_attribution.py
echo "no-motif attribution done -- trends in no_motif_attr_output/"
