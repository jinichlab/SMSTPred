#!/bin/bash
#SBATCH --job-name=find_bound_struct
#SBATCH --ntasks=1
#SBATCH --output=find_bound_structures_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00

# Due-diligence for the apo/AlphaFold-only motifs: for each, expand the FIMO-matching
# reps to their full 50%-id clusters and ask whether any cluster MEMBER carrying the
# motif has a substrate/cofactor-bound PDB we could use instead of the apo/AF model.
#
# NOTE: this job makes outbound HTTPS calls (rest.uniprot.org id-mapping + data.rcsb.org).
# It must run on a node/partition with internet access; if compute nodes are firewalled,
# run it on the login node instead:  bash run_find_bound_structures.sh
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

python find_bound_structures.py
