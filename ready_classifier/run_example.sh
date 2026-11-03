#!/bin/bash
#SBATCH --job-name=smstpred_example
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=run_example_%j.log
# Run the classifier on the six example sequences (a GPU is optional: remove the --gpus line and add `--device cpu` to run without one).
#   conda env create -f environment.yml && conda activate smstpred-classify     (once)
#   bash run_example.sh        or      sbatch run_example.sh
cd "$(dirname "$0")" || exit 1
python classify.py example/example.fasta -o example/example_predictions.csv
