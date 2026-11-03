#!/bin/bash
#SBATCH --job-name=classical_ml
#SBATCH --ntasks=1
#SBATCH --array=1-16
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=run_classical_%A_%a.log

# Stage 2: the classical-ML arm of the comparison (Fig S4): XGBoost, Random Forest, Logistic Regression and four SVMs on the ProtT5 embeddings, on the 16 splits of the bake-off, one split per array
# task (classical_ml.py). CPU only; the SVMs (probability=True) dominate the time. Appends iterNN_classical_* to <release>/benchmark/iterations/; run_bakeoff_report.sh then includes them.
#   sbatch run_classical.sh                    # all 16 splits
#   sbatch --array=1 --export=ALL,SMOKE=1 run_classical.sh      # tiny end-to-end test into <release>/benchmark/_smoke
#   CLASSICAL_METHODS="XGBoost,Random Forest" sbatch run_classical.sh      # only some classifiers
# Release folder: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR. Needs the splits (make splits) and the embeddings of the release. Resumable.
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
python -u classical_ml.py --iteration "${SLURM_ARRAY_TASK_ID:-${ITERATION:?run with sbatch (SLURM_ARRAY_TASK_ID) or set ITERATION}}" ${SMOKE:+--smoke}
