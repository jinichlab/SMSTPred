#!/bin/bash
#SBATCH --job-name=mmseqs_clusters
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=run_mmseqs_clusters_output.log

# Stage 1: MMseqs2 clusters (30 % identity) and cluster-disjoint splits of the typed sets (SAM benchmark: also the 16 train / validation / test splits of the bake-off, tvt_splits.py). CPU only, no network. Run after 0_uniprot_pull/run_bin.sh.
#   sbatch run_mmseqs_clusters.sh                                  # the two SAM clusterings (benchmark <= 1,000 aa, deployment <= 2,000 aa): about a minute
#   CLUSTER_WHAT=up sbatch --cpus-per-task=32 --mem=100G run_mmseqs_clusters.sh      # UP (+ the SAM benchmark proteins), ~210,000 sequences: the big one
#   CLUSTER_WHAT=sam-benchmark | sam-deployment | up | sam (both SAM ones, default)
#   CLUSTER_ARGS="--min-seq-id 0.4 --tag-suffix _x" sbatch run_mmseqs_clusters.sh    # another identity; the tag keeps the variants apart
# Reads $UNIPROT_PULL_DIR (default <repo>/data/uniprot_pull)/release_*/ (labels_SAM.csv, labels_UP.csv with their in_* flags, and the pulled tables for the sequences);
# writes <release>/clusters/ and the cluster_* columns of the label tables. mmseqs must be on the PATH (conda: bioconda mmseqs2) or in $MMSEQS.
# Environment: the conda environment that is already active, else the one `make env` builds (SMSTPRED_ENV renames it); it needs mmseqs, numpy, pandas.
set -e

here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
case "${CLUSTER_WHAT:-sam}" in
    sam)            python -u mmseqs_clusters.py sam --subset benchmark ${CLUSTER_ARGS:-}; python -u tvt_splits.py; python -u mmseqs_clusters.py sam --subset deployment ${CLUSTER_ARGS:-} ;;
    sam-benchmark)  python -u mmseqs_clusters.py sam --subset benchmark ${CLUSTER_ARGS:-}; python -u tvt_splits.py ;;
    sam-deployment) python -u mmseqs_clusters.py sam --subset deployment ${CLUSTER_ARGS:-} ;;
    up)             python -u mmseqs_clusters.py up ${CLUSTER_ARGS:-} ;;
    *) echo "CLUSTER_WHAT must be sam | sam-benchmark | sam-deployment | up" >&2; exit 1 ;;
esac
