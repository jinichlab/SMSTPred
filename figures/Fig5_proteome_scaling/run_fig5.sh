#!/bin/bash
#SBATCH --job-name=fig5
#SBATCH --ntasks=1
#SBATCH --output=run_fig5_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00

# Fig 5 (SAM-MTase count vs proteome size) for a refresh release: build the per-organism table, then the figure + fits.
# Needs <release>/predictions_all_proteomes.csv.gz (python 3_predictions/combine_predictions.py).   RELEASE_DIR=<release> sbatch figures/Fig5_proteome_scaling/run_fig5.sh   Output: <release>/fig5/
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
python build_proteome_vs_sam_table.py
REL="$(python -c "import os,sys;sys.path.insert(0,'../../1_clustering');from pathlib import Path;from mmseqs_clusters import newest_release;print(os.environ.get('RELEASE_DIR') or newest_release(Path(os.environ.get('UNIPROT_PULL_DIR','../../data/uniprot_pull')),os.environ.get('RELEASE')))")"
python proteome_vs_sam_domain_highlight.py --dir "$REL/fig5"
python scaling_loglog.py --dir "$REL/fig5"                       # Fig 5 as power-law fits (exponents, level differences, clade-bootstrap CIs): fig5_scaling_loglog.png, fig5_scaling_numbers.csv
python ../../4_phylo_trends/lineage_statistics.py --release "$REL"      # the lineage comparisons quoted in the Results (<release>/fig5/lineage_statistics.csv/.txt)
