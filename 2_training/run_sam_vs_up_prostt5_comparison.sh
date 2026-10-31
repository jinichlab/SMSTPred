#!/bin/bash
#SBATCH --job-name=sam_vs_up
#SBATCH --ntasks=1
#SBATCH --output=sam_vs_up_prostt5_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=06:00:00

set -e

# SAM-trained vs UP-trained ProstT5 MLP on a held-out SAM test set, over N=16
# repeated cluster-disjoint 90/10 splits (matches Panel A's design). MVP: only
# the exact test entries are removed from UP training (homolog exclusion deferred). The UP model trains for 400 epochs
# with the per-label losses set in UP_CFG (sam_vs_up_prostt5_comparison.py).
# Inputs are read from <repo>/data, or from $SMSTPRED_DATA.
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
python -u sam_vs_up_prostt5_comparison.py

# Stats + standalone Panel B figure (CPU only).
python sam_vs_up_prostt5_stats.py     # Nadeau-Bengio corrected paired test (rho=1/9)
python plot_sam_vs_up_prostt5.py      # standalone Panel B with significance brackets
python plot_fig2_main.py              # refresh the Fig 2 composite (A over B)
