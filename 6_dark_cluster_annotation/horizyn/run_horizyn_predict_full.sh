#!/bin/bash
#SBATCH --job-name=horizyn_full
#SBATCH --gres=gpu:rtx6000ada:1
#SBATCH --output=horizyn_predict_full.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00

# Score ALL 3,478 dark-cluster proteins against ALL 36,552 Rhea reactions
# (every EC class) with the pretrained Horizyn INF checkpoint. No training.
# Est. ~15 min reaction fingerprinting (CPU) + ~30-60 min protein embedding
# (GPU, ProtT5-XL) + seconds for the score matmul.
# The python scripts of this folder import the `horizyn` package from the folder they sit in: copy them (here run_horizyn_predict_full.py and embed_prott5.py) into $HORIZYN_DIR first.
# The Horizyn environment is not part of environment.yml: it comes with its own install (the versions used here are in horizyn/environment_horizyn.yml).
# Uses the active environment if it is horizyn, else activates ${HORIZYN_ENV:-horizyn}.
if [ "${CONDA_DEFAULT_ENV:-}" != "${HORIZYN_ENV:-horizyn}" ] && command -v conda > /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "${HORIZYN_ENV:-horizyn}" || { echo "conda environment ${HORIZYN_ENV:-horizyn} not found: create it from horizyn/environment_horizyn.yml (see 6_dark_cluster_annotation/README.md)" >&2; exit 1; }
fi
cd "${HORIZYN_DIR:?set HORIZYN_DIR to your Horizyn checkout (the folder holding checkpoints/ and configs/)}"
python run_horizyn_predict_full.py
