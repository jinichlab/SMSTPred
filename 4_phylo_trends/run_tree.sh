#!/bin/bash
#SBATCH --job-name=tree
#SBATCH --ntasks=1
#SBATCH --output=run_tree_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

# Fig 3 for a refresh release, with no iTOL session: combine the per-proteome predictions, build the taxonomy tree + annotation files, draw the figure.
#   RELEASE_DIR=<release> sbatch 4_phylo_trends/run_tree.sh        (default release: the newest in $UNIPROT_PULL_DIR)
# Needs the proteome predictions of Stage 3 (<release>/predictions/*.csv). Output: <release>/tree/ (iTOL files, tree_stats.csv, fig3_tree.png/.svg/.pdf).
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
python ../3_predictions/combine_predictions.py
python reference_tree.py --release
REL="$(python -c "import os,sys;sys.path.insert(0,'../1_clustering');from pathlib import Path;from mmseqs_clusters import newest_release;print(os.environ.get('RELEASE_DIR') or newest_release(Path(os.environ.get('UNIPROT_PULL_DIR','../data/uniprot_pull')),os.environ.get('RELEASE')))")"
python ../figures/Fig3_tree_landscape/draw_tree.py --src "$REL/tree"
