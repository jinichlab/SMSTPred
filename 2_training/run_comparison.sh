#!/bin/bash
#SBATCH --job-name=sam_clf_cmp
#SBATCH --ntasks=1
#SBATCH --output=run_comparison_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G

set -e
source ~/miniconda3/etc/profile.d/conda.sh

# Sweep: 16 iters x {cluster, random} x {full, balanced} x 16 classifiers.
conda activate classifier_env
python sam_full_vs_balanced_classifier_comparison.py

# Stats + figures off the two CSVs above. scikit_posthocs lives in base, not
# classifier_env, so the downstream steps run there (CPU only, no GPU needed).
conda activate base
python classifier_mcc_stats.py      # Friedman + Bonferroni-Dunn vs MLP_ProstT5
python plot_cluster_mcc_zoom.py     # Fig 2A standalone + balanced supplementary
python build_metrics_table.py       # Table 2: 16 classifiers x 5 metrics
python plot_fig2_main.py            # Fig 2 composite (A over B)
