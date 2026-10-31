#!/bin/bash
#SBATCH --job-name=bakeoff_mlp
#SBATCH --ntasks=1
#SBATCH --output=run_bakeoff_mlp_only_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00

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
# Re-train only the six MLP variants of the classifier comparison (full training set, 300 epochs, the loss map in
# MLP_VARIANTS; about 3 minutes on an RTX 6000 Ada). The train/test splits are identical to a full run of
# run_comparison.sh, because the Profile-HMM's random draws are replayed, so these rows can be dropped into the
# full results with merge_bakeoff_rows.py. Inputs are read from <repo>/data, or from $SMSTPRED_DATA.
export BAKEOFF_ONLY=mlp
export BAKEOFF_OUT=sam_classifier_comparison_output_mlp_only

mkdir -p "$BAKEOFF_OUT"
cp -n sam_classifier_comparison_output/full_clusters_t0.6.pkl "$BAKEOFF_OUT"/     # cached clusters: same splits as the committed run

python -u sam_classifier_comparison.py
