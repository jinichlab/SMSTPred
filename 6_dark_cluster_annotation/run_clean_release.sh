#!/bin/bash
#SBATCH --job-name=clean_dark
#SBATCH --gpus=1
#SBATCH --output=clean_dark.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# CLEAN (maxsep + p-value) on the release's dark proteins; no arguments. Needs <release>/dark_annotation/<name>.fasta (run_dark_inputs.sh).
#   sbatch 6_dark_cluster_annotation/run_clean_release.sh
# The CLEAN checkout is annotation.clean_app_dir of refresh_config.json (or $CLEAN_APP_DIR); its conda environment annotation.clean_env (or $CLEAN_ENV).
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
source "$here/../workflow/stage_env.sh"
set -e
A="$REL/dark_annotation"
export DARK_NAME="$(cat "$A/dark_name.txt")"
export CLEAN_APP_DIR="$(cfg clean_app_dir)"; [ -n "$CLEAN_APP_DIR" ] || { echo "set annotation.clean_app_dir in refresh_config.json (or CLEAN_APP_DIR)" >&2; exit 1; }
export CLEAN_ENV="$(cfg clean_env)"; export CLEAN_ENV="${CLEAN_ENV:-clean}"
export DARK_FASTA="$A/$DARK_NAME.fasta" DARK_OUT="$A/clean"
conda deactivate 2>/dev/null || true       # run_clean_dark.sh activates the CLEAN environment itself
bash "$REPO/6_dark_cluster_annotation/clean/run_clean_dark.sh"
