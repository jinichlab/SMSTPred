#!/bin/bash
#SBATCH --job-name=motif_dedup
#SBATCH --ntasks=1
#SBATCH --output=dedup_output.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Cluster the 106k single-label SAM-MTase sequences at 50% identity so motif
# discovery sees one representative per family (not thousands of paralogs).
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

rm -rf tmp_mmseqs
mmseqs easy-cluster corpus_all.fasta clust tmp_mmseqs \
    --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads 16
python 2_make_reps.py
