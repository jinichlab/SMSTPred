#!/bin/bash
#SBATCH --job-name=clean_dark
#SBATCH --gpus=1
#SBATCH --output=clean_dark.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# CLEAN (maxsep, then p-value; the second reuses the ESM embeddings of the first) for the dark proteins of a named FASTA.
#   DARK_NAME=dark53 CLEAN_APP_DIR=<CLEAN checkout>/app DARK_FASTA=<release>/dark_annotation/dark53.fasta DARK_OUT=<release>/dark_annotation/clean sbatch run_clean_dark.sh
# DARK_NAME (default dark33) is the file stem in data/inputs/ of the CLEAN checkout; DARK_FASTA is copied there; the two result files are copied to DARK_OUT.
# The CLEAN environment is not part of environment.yml (versions in clean/environment_clean.yml); the active one is used if it is clean, else ${CLEAN_ENV:-clean}.
if [ "${CONDA_DEFAULT_ENV:-}" != "${CLEAN_ENV:-clean}" ] && command -v conda > /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "${CLEAN_ENV:-clean}" || { echo "conda environment ${CLEAN_ENV:-clean} not found: create it from clean/environment_clean.yml" >&2; exit 1; }
fi
HERE="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && HERE="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
HERE="$(cd "$HERE" && pwd)"
NAME="${DARK_NAME:-dark33}"
cd "${CLEAN_APP_DIR:?set CLEAN_APP_DIR to the app/ directory of your CLEAN checkout}"
export PYTHONPATH="$CLEAN_APP_DIR/esm:$PYTHONPATH"      # the vendored fair-esm (the env's site-packages esm is ESM-3, wrong API)
[ -n "${DARK_FASTA:-}" ] && cp "$DARK_FASTA" "data/inputs/$NAME.fasta"
set -e
python -u CLEAN_infer_fasta.py -d "$NAME"
cp "$HERE/CLEAN_infer_fasta_pvalue.py" ./CLEAN_infer_fasta_pvalue_named.py
python -u CLEAN_infer_fasta_pvalue_named.py -d "$NAME"
if [ -n "${DARK_OUT:-}" ]; then
    mkdir -p "$DARK_OUT"
    cp "results/inputs/${NAME}_maxsep.csv" "results/inputs/${NAME}_pvalue.csv" "$DARK_OUT/"
fi
