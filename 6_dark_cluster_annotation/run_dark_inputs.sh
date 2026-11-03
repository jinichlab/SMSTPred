#!/bin/bash
#SBATCH --job-name=dark_inputs
#SBATCH --ntasks=1
#SBATCH --output=run_dark_inputs_output.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# The dark proteins of the release's map (5_umap/run_clusters.sh) as CLEAN / Horizyn inputs: <release>/dark_annotation/<name>.fasta + <name>_members.csv + dark_name.txt.
#   sbatch 6_dark_cluster_annotation/run_dark_inputs.sh
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
source "$here/../workflow/stage_env.sh"
set -e
cd "$here"
python make_dark_inputs.py --release "$REL"
