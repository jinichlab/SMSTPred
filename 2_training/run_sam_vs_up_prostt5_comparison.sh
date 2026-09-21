#!/bin/bash
#SBATCH --job-name=sam_vs_up
#SBATCH --ntasks=1
#SBATCH --output=sam_vs_up_prostt5_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=06:00:00

set -e

# SAM-trained vs UP-trained ProstT5 MLP on a held-out SAM test set, over N=16
# repeated cluster-disjoint 90/10 splits (matches Panel A's design). MVP: only
# the exact test entries are removed from UP training (homolog exclusion deferred).
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env
python -u sam_vs_up_prostt5_comparison.py

# Stats + standalone Panel B figure (CPU only). scipy/matplotlib live in base.
conda activate base
python sam_vs_up_prostt5_stats.py     # Nadeau-Bengio corrected paired test (rho=1/9)
python plot_sam_vs_up_prostt5.py      # standalone Panel B with significance brackets
python plot_fig2_main.py              # refresh the Fig 2 composite (A over B)
