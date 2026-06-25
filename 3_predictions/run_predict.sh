#!/bin/bash
#SBATCH --job-name=sam_predict
#SBATCH --ntasks=1
#SBATCH --output=predict_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Usage: edit --input to point at your ProstT5 embedding H5 file or directory.
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python predict.py --input "$1" --out-dir predictions_output
