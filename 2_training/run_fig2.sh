#!/bin/bash
#SBATCH --job-name=fig2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=run_fig2_%j.log
# Fig 2 (A and B in one figure) of a release: after run_bakeoff_report.sh (A, Table 1) and run_sam_vs_up_report.sh (B). Writes <release>/benchmark/Fig2.png.
#   sbatch run_fig2.sh        (or: bash run_fig2.sh)      Release folder: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR.
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
python -u plot_fig2_release.py
