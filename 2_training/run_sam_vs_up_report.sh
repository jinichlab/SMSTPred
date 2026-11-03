#!/bin/bash
#SBATCH --job-name=sam_vs_up_report
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=run_sam_vs_up_report_%j.log

# Stage 2: after run_sam_vs_up.sh -- merge the 16 splits, the paired tests (Wilcoxon, paired t, Nadeau-Bengio) of SAM-trained vs UP-trained, the figure and the training curves of the UP arm
# (sam_vs_up_report.py). CPU, seconds.   sbatch run_sam_vs_up_report.sh   (or bash run_sam_vs_up_report.sh)
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
python -u sam_vs_up_report.py
