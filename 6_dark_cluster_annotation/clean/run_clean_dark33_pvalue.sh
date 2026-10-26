#!/bin/bash
#SBATCH --job-name=clean_dark33_pval
#SBATCH --gpus=1
#SBATCH --output=clean_dark33_pvalue.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate clean
cd "${CLEAN_APP_DIR:?set CLEAN_APP_DIR to the app/ directory of your CLEAN checkout}"
export PYTHONPATH="$CLEAN_APP_DIR/esm:$PYTHONPATH"
python -u CLEAN_infer_fasta_pvalue.py
