#!/bin/bash
#SBATCH --job-name=bakeoff_report
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=run_bakeoff_report_%j.log

# Stage 2: after run_bakeoff.sh -- merge the 16 splits, the statistics (Friedman, Nemenyi, Bonferroni-Dunn, Wilcoxon), Fig 2A, Table 1, the paired tests against the search baselines
# and the train / validation loss figures with the epochs early stopping kept (bakeoff_report.py). CPU, about a minute.
#   sbatch run_bakeoff_report.sh        (or: bash run_bakeoff_report.sh)
# Release folder: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR. Writes into <release>/benchmark/ (BAKEOFF_OUT overrides).
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
python -u bakeoff_report.py
