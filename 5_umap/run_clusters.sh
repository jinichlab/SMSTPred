#!/bin/bash
#SBATCH --job-name=clusters
#SBATCH --ntasks=1
#SBATCH --output=run_clusters_output.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Dark / metazoan-protein / vertebrate-DNA clusters of a release's UMAP, and the Fig 7 grids. Needs <release>/umap/ from global_umap.py.
#   RELEASE_DIR=<release> sbatch 5_umap/run_clusters.sh      Output: <release>/umap/ (cluster tables, crosswalk, report) and <release>/umap/fig7/
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
set -e
python test_cluster_selection.py
python cluster_selection.py
python ../figures/Fig7_clusters/build_cluster_grids.py --release
python ../figures/Fig7_clusters/make_fig7_clusters.py --dir "$(python ../workflow/release_dir.py)/umap/fig7"        # the composite (run_tables.sh of Stage 6 draws it again with the Fig 8 clusters outlined)
