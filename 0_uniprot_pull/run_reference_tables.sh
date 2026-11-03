#!/bin/bash
#SBATCH --job-name=ref_tables
#SBATCH --ntasks=1
#SBATCH --output=reference_tables_output.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
# The public reference tables of the dark-cluster annotation (ExPASy ENZYME, Rhea rhea2ec + rhea-directions) into <release>/reference_tables/ (network).
#   sbatch 0_uniprot_pull/run_reference_tables.sh      (release: $RELEASE_DIR, else the newest in $UNIPROT_PULL_DIR; a file already there is kept)
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
source "$here/../workflow/stage_env.sh"
set -e
cd "$here"
python fetch_reference_tables.py --release "$(basename "$REL" | sed 's/^release_//')"
