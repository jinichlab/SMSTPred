#!/bin/bash
#SBATCH --job-name=umap
#SBATCH --ntasks=1
#SBATCH --output=run_umap_output.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Global UMAP of a refresh release from its ProtT5 embeddings (115k proteins): coordinates, HDBSCAN clusters, static + interactive map, neighbour purity.
# Needs the proteome embeddings (embeddings/run_embeddings.sh proteomes) and <release>/predictions_all_proteomes.csv.gz (3_predictions/combine_predictions.py).
#   RELEASE_DIR=<release> sbatch 5_umap/run_umap.sh      Output: <release>/umap/
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
python global_umap.py
