#!/bin/bash
#SBATCH --job-name=fig8
#SBATCH --ntasks=1
#SBATCH --output=run_fig8_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
# Refreshed Fig 8 from a release: RELEASE_DIR=<release> [SMSTPRED_ENV=smstpred_fresh] sbatch run_fig8.sh
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${SMSTPRED_ENV:-smstpred_fresh}"
python make_fig8_reaction_profiles.py --release "${RELEASE_DIR:?set RELEASE_DIR}"
