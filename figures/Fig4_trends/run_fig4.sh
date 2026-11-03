#!/bin/bash
#SBATCH --job-name=fig4
#SBATCH --ntasks=1
#SBATCH --output=run_fig4_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Fig 4 (substrate-load trends) for a refresh release: the four violin panels, the 2x2 composite, fig4_statistics.csv. Needs <release>/predictions_all_proteomes.csv.gz
# (python 3_predictions/combine_predictions.py).   RELEASE_DIR=<release> sbatch figures/Fig4_trends/run_fig4.sh      Output: <release>/fig4/
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
python build_fig4_panels.py --release
REL="$(python -c "import os,sys;sys.path.insert(0,'../../1_clustering');from pathlib import Path;from mmseqs_clusters import newest_release;print(os.environ.get('RELEASE_DIR') or newest_release(Path(os.environ.get('UNIPROT_PULL_DIR','../../data/uniprot_pull')),os.environ.get('RELEASE')))")"
python make_fig4_trends.py --dir "$REL/fig4"
