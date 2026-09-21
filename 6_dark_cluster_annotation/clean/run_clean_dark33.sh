#!/bin/bash
#SBATCH --job-name=clean_dark33
#SBATCH --gpus=1
#SBATCH --output=clean_dark33.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate clean
cd /home/nroethler/Code/abhiram/CLEAN/app
# use the vendored fair-esm (env's site-packages esm is EvolutionaryScale ESM-3, wrong API)
export PYTHONPATH=/home/nroethler/Code/abhiram/CLEAN/app/esm:$PYTHONPATH
python -u CLEAN_infer_fasta.py -d dark33
