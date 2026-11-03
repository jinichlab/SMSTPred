#!/bin/bash
#SBATCH --job-name=domain_table
#SBATCH --ntasks=1
#SBATCH --output=domain_table_output.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Table S3 for a refresh release: corpus (labelled SAM enzymes + the ProtT5-ensemble proteome calls) -> InterPro enrichment -> table.
#   RELEASE_DIR=<release> sbatch run_domain_table.sh      (default release: the newest in $UNIPROT_PULL_DIR)
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
set -e
python 1_build_corpus.py --release
REL="$(python -c "import os,sys;sys.path.insert(0,'../1_clustering');from pathlib import Path;from mmseqs_clusters import newest_release;print(os.environ.get('RELEASE_DIR') or newest_release(Path(os.environ.get('UNIPROT_PULL_DIR','../data/uniprot_pull')),os.environ.get('RELEASE')))")"
python 3_domain_enrichment.py --dir "$REL/motif_analysis"
python ../figures/TableS3_domain_enrichment/build_domain_enrichment_table.py --src "$REL/motif_analysis/domain_enrichment.csv" --out "$REL/motif_analysis"
