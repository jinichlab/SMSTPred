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
source ~/miniconda3/etc/profile.d/conda.sh
conda activate horizyn
cd "${HORIZYN_DIR:?set HORIZYN_DIR to your Horizyn checkout (the folder holding checkpoints/ and configs/)}"
python run_horizyn_predict_full.py
