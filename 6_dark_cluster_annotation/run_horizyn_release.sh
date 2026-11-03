#!/bin/bash
#SBATCH --job-name=horizyn_dark
#SBATCH --gres=gpu:rtx6000ada:1
#SBATCH --output=horizyn_dark.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
# Horizyn (all 36,552 Rhea reactions) on the release's dark proteins; no arguments. Needs <release>/dark_annotation/<name>.fasta + <name>_members.csv (run_dark_inputs.sh).
#   sbatch 6_dark_cluster_annotation/run_horizyn_release.sh
# The Horizyn checkout is annotation.horizyn_dir of refresh_config.json (or $HORIZYN_DIR; the folder holding checkpoints/, configs/ and horizyn/), its environment annotation.horizyn_env.
# The checkout is never written to: the run folder <release>/dark_annotation/horizyn_run/ links its pieces and holds copies of the two scripts and all outputs.
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
source "$here/../workflow/stage_env.sh"
set -e
A="$REL/dark_annotation"; RUN="$A/horizyn_run"
NAME="$(cat "$A/dark_name.txt")"
CK="$(cfg horizyn_dir)"; [ -d "$CK" ] || { echo "set annotation.horizyn_dir in refresh_config.json (or HORIZYN_DIR) to your Horizyn checkout" >&2; exit 1; }
mkdir -p "$RUN"
for x in checkpoints configs horizyn; do ln -sfn "$CK/$x" "$RUN/$x"; done
ln -sfn "$REPO/6_dark_cluster_annotation/horizyn/inputs" "$RUN/inputs"
cp "$REPO/6_dark_cluster_annotation/horizyn/run_horizyn_predict_full.py" "$REPO/6_dark_cluster_annotation/horizyn/embed_prott5.py" "$RUN/"
export HORIZYN_DIR="$RUN" SMSTPRED_DARK_FASTA="$A/$NAME.fasta" SMSTPRED_DARK_MEMBERS="$A/${NAME}_members.csv" SMSTPRED_STAGE6="$REPO/6_dark_cluster_annotation"
export HORIZYN_ENV="$(cfg horizyn_env)"; export HORIZYN_ENV="${HORIZYN_ENV:-horizyn}"
conda deactivate 2>/dev/null || true
bash "$REPO/6_dark_cluster_annotation/horizyn/run_horizyn_predict_full.sh" 2>&1 | tee "$RUN/horizyn_predict_full.log"
test -s "$RUN/horizyn_full_topk.csv"      # the run wrote its result (a failure above leaves no file)
