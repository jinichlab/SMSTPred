#!/bin/bash
#SBATCH --job-name=clean_dark33
#SBATCH --gpus=1
#SBATCH --output=clean_dark33.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# The CLEAN environment is not part of environment.yml: it comes with its own install (the versions used here are in clean/environment_clean.yml).
# Uses the active environment if it is clean, else activates ${CLEAN_ENV:-clean}.
if [ "${CONDA_DEFAULT_ENV:-}" != "${CLEAN_ENV:-clean}" ] && command -v conda > /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "${CLEAN_ENV:-clean}" || { echo "conda environment ${CLEAN_ENV:-clean} not found: create it from clean/environment_clean.yml (see 6_dark_cluster_annotation/README.md)" >&2; exit 1; }
fi
cd "${CLEAN_APP_DIR:?set CLEAN_APP_DIR to the app/ directory of your CLEAN checkout}"
# use the vendored fair-esm (env's site-packages esm is EvolutionaryScale ESM-3, wrong API)
export PYTHONPATH="$CLEAN_APP_DIR/esm:$PYTHONPATH"
python -u CLEAN_infer_fasta.py -d dark33
