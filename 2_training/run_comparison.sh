#!/bin/bash
#SBATCH --job-name=sam_clf_cmp
#SBATCH --ntasks=1
#SBATCH --output=run_comparison_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G

set -e

# Environment: the conda environment that is already active, else the one `make env` builds (SMSTPRED_ENV renames it).
# Runs as `sbatch run_*.sh` (the job inherits the submitting shell's environment) or as plain `bash run_*.sh`.
# The scripts use relative paths, so start in this script's folder (under sbatch $0 is a spooled copy: ask SLURM where it came from).
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
# Sweep: 16 iters x {cluster, random} x 16 classifiers, trained on the full SAM set (add BAKEOFF_DATASETS=full,balanced for the optional balanced arm).
python -u sam_classifier_comparison.py

# Stats + figures off the two CSVs above (CPU only, no GPU needed).
python classifier_mcc_stats.py      # Friedman + Bonferroni-Dunn vs MLP_ProstT5
python plot_cluster_mcc_zoom.py     # Fig 2A and Fig S4
python build_metrics_table.py       # Table 1: nine representation and search-tool methods x 5 metrics
python plot_fig2_main.py            # Fig 2 composite (A over B)
