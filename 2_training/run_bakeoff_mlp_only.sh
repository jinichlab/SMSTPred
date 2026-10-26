#!/bin/bash
#SBATCH --job-name=bakeoff_mlp
#SBATCH --ntasks=1
#SBATCH --output=run_bakeoff_mlp_only_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00

set -e

# Re-train only the six MLP variants of the classifier comparison (full training set, 300 epochs, the loss map in
# MLP_VARIANTS; about 15 minutes on one GPU). The train/test splits are identical to a full run of
# run_comparison.sh, because the Profile-HMM's random draws are replayed, so these rows can be dropped into the
# full results with merge_bakeoff_rows.py. Inputs are read from <repo>/data, or from $SMSTPRED_DATA.
export BAKEOFF_ONLY=mlp
export BAKEOFF_OUT=sam_classifier_comparison_output_mlp_only

mkdir -p "$BAKEOFF_OUT"
cp -n sam_classifier_comparison_output/full_clusters_t0.6.pkl "$BAKEOFF_OUT"/     # cached clusters: same splits as the committed run

source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env
python -u sam_classifier_comparison.py
