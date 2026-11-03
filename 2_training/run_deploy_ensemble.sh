#!/bin/bash
#SBATCH --job-name=deploy_ensemble
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=run_deploy_ensemble_%j.log

# Stage 2: the deployment model of the refresh: a 6-member ProtT5 MLP ensemble from a 6-fold CLUSTER partition of the deployment set (deploy_ensemble.py): about 5 minutes on one GPU.
#   sbatch run_deploy_ensemble.sh                       # writes <release>/deployment/ (the files 3_predictions/predict.py reads: DEPLOY_DIR=<release>/deployment)
#   SMOKE=1 sbatch run_deploy_ensemble.sh               # tiny end-to-end test into <release>/deployment/_smoke
# Needs: labels_SAM.csv with the deployment clusters (1_clustering/run_mmseqs_clusters.sh) and embeddings/embeddings_SAM_ProtT5.h5 (embeddings/run_embeddings.sh).
# Release folder: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR. Settings: ../refresh_config.json ("deployment").
set -e

here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
python -u test_deploy_ensemble.py
python -u deploy_ensemble.py ${SMOKE:+--smoke}
