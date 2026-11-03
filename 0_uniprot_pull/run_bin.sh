#!/bin/bash
#SBATCH --job-name=bin_substrates
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=run_bin_output.log

# Stage 0: type the pulled proteins by substrate (Metabolite / Protein / RNA / DNA) from their reaction annotation. CPU only, no network.
#   sbatch run_bin.sh                                  # sam + up of the newest pull (labels_SAM.csv, labels_UP.csv, report, rule list)
#   BIN_WHAT=sam sbatch run_bin.sh                     # one part: sam | up | all
#   BIN_ARGS="--deployment-max-len 1500" sbatch run_bin.sh
#   BIN_ARGS="--benchmark-sources reaction_curated reaction_uncurated" sbatch run_bin.sh   # which label sources count as in_benchmark
# Reads $UNIPROT_PULL_DIR (default <repo>/data/uniprot_pull)/release_*/ and the archived labels from $SMSTPRED_DATA (default <repo>/data).
# Environment: the conda environment that is already active, else the one `make env` builds (SMSTPRED_ENV renames it).
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
python -u test_substrate_binning.py
# BIN_ARGS="--refresh-rhea" downloads a new Rhea compound table first (network; the repo carries the 2026-11-02 one)
python -u bin_substrates.py "${BIN_WHAT:-all}" ${BIN_ARGS:-}
