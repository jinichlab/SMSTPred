#!/bin/bash
#SBATCH --job-name=motif_dedup
#SBATCH --ntasks=1
#SBATCH --output=dedup_output.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Cluster the 106k single-label SAM-MTase sequences at 50% identity so motif
# discovery sees one representative per family (not thousands of paralogs).
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

rm -rf tmp_mmseqs
mmseqs easy-cluster corpus_all.fasta clust tmp_mmseqs \
    --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads 16
python 2_make_reps.py
