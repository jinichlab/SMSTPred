#!/bin/bash
#SBATCH --job-name=uniprot_pull
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=run_pull_output.log

# Stage 0: pull the current UniProt data through the UniProt REST API (network access needed; compute nodes with internet are fine).
#   sbatch run_pull.sh                              # sam + up + the 1,259 proteomes, then the comparison with the archived data
#   PULL_WHAT=sam sbatch run_pull.sh                # only one part: sam | up | proteomes | all
#   PULL_WHAT=resplit PULL_ARGS="--length-view 1500" sbatch run_pull.sh   # change the length cutoff of the newest pull (no download)
#   PULL_ARGS="--limit 3" sbatch run_pull.sh        # smoke test (3 proteomes, a few hundred rows; written under <out>/_smoke/)
#   PULL_ARGS="--expect-release 2026_03" ...        # stop unless UniProt serves that release
# Output goes to $UNIPROT_PULL_DIR (default <repo>/data/uniprot_pull)/release_<UniProt release>/ and never overwrites archived data.
# The proteome part is resumable (proteome_status.jsonl): submit the same job again after an interruption.
# The comparison reads the archived tables from $SMSTPRED_DATA (default <repo>/data).
# Environment: the conda environment that is already active, else the one `make env` builds (SMSTPRED_ENV renames it).
# Runs as `sbatch run_*.sh` (the job inherits the submitting shell's environment) or as plain `bash run_*.sh`.
# The scripts use relative paths, so start in this script's folder (under sbatch $0 is a spooled copy: ask SLURM where it came from).
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
python -u pull_uniprot.py "${PULL_WHAT:-all}" ${PULL_ARGS:-}

# what changed against the archived data (sets, sequences, proteome sizes); only after a complete pull (not after a --limit smoke test)
if { [ "${PULL_WHAT:-all}" = all ] || [ "${PULL_WHAT:-all}" = resplit ]; } && [[ " ${PULL_ARGS:-} " != *" --limit "* ]]; then
    python -u compare_to_previous.py ${COMPARE_ARGS:-}
fi
