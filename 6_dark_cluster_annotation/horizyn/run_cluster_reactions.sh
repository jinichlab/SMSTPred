#!/bin/bash
#SBATCH --job-name=horizyn_rxn_cluster
#SBATCH --gres=gpu:rtx6000ada:1
#SBATCH --output=cluster_reactions.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00

# Cluster all 36,552 Rhea reactions in Horizyn's own query-embedding space
# (HDBSCAN), so protein clusters can be labeled by dominant REACTION FAMILY
# rather than a single raw Rhea ID (which can artificially fragment one real
# chemistry, as with c758's menaquinone chain-length variants).
# The python scripts of this folder import the `horizyn` package from the folder they sit in: copy them (here cluster_rhea_reactions.py) into $HORIZYN_DIR first.
# The Horizyn environment is not part of environment.yml: it comes with its own install (the versions used here are in horizyn/environment_horizyn.yml).
# Uses the active environment if it is horizyn, else activates ${HORIZYN_ENV:-horizyn}.
if [ "${CONDA_DEFAULT_ENV:-}" != "${HORIZYN_ENV:-horizyn}" ] && command -v conda > /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "${HORIZYN_ENV:-horizyn}" || { echo "conda environment ${HORIZYN_ENV:-horizyn} not found: create it from horizyn/environment_horizyn.yml (see 6_dark_cluster_annotation/README.md)" >&2; exit 1; }
fi
cd "${HORIZYN_DIR:?set HORIZYN_DIR to your Horizyn checkout (the folder holding checkpoints/ and configs/)}"
python cluster_rhea_reactions.py
