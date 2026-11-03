#!/bin/bash
#SBATCH --job-name=sam_predict_rel
#SBATCH --ntasks=1
#SBATCH --output=predict_release_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
# Stage 3 for a refresh release, no arguments: score the release's reference-proteome embeddings (ProtT5) with its deployed ensemble, then combine the per-proteome files.
#   sbatch 3_predictions/run_predict_release.sh        (release: $RELEASE_DIR, else the newest in $UNIPROT_PULL_DIR)
# Needs <release>/deployment/ (2_training/run_deploy_ensemble.sh) and <release>/embeddings/proteomes_ProtT5/ (embeddings/run_embeddings.sh with EMBED_WHAT=proteomes).
# Writes <release>/predictions/<proteome>.csv and <release>/predictions_all_proteomes.csv.gz (+ summary).
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
source "$here/../workflow/stage_env.sh"
set -e
cd "$here"
for need in "$REL/deployment" "$REL/embeddings/proteomes_ProtT5"; do [ -e "$need" ] || { echo "ERROR: $need not found: run the earlier stages first" >&2; exit 1; }; done
DEPLOY_DIR="$REL/deployment" python predict.py --input "$REL/embeddings/proteomes_ProtT5" --out-dir "$REL/predictions"
python combine_predictions.py
