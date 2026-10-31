#!/bin/bash
#SBATCH --job-name=verify_upgrade
#SBATCH --ntasks=1
#SBATCH --output=verify_upgrade_motifs_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# Confirm each apo/AF-motif upgrade structure's protein actually carries the target
# motif: fetch its UniProt sequence and FIMO-scan with the class's STREME motifs.
# Needs outbound HTTPS (rest.uniprot.org) + fimo (MEME suite, in the smstpred environment).
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

python verify_upgrade_motifs.py
