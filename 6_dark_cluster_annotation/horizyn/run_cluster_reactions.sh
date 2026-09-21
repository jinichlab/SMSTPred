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
source ~/miniconda3/etc/profile.d/conda.sh
conda activate horizyn
cd /home/nroethler/Code/horizyn/horizyn_local
python cluster_rhea_reactions.py
