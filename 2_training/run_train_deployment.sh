#!/bin/bash
#SBATCH --job-name=deploy_train
#SBATCH --ntasks=1
#SBATCH --output=train_deployment_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python train_deployment_ensemble.py
