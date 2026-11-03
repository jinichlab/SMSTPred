#!/bin/bash
#SBATCH --job-name=no_motif_umap
#SBATCH --ntasks=1
#SBATCH --output=no_motif_umap_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
# Fig 11 of a release: the UMAP regions without a motif. Needs <release>/motif_analysis/motif_status_all.csv (run_fimo_all.sh) and <release>/umap/ (5_umap/run_umap.sh).
#   [MOTIF_DIR=<release>/motif_analysis | auto] sbatch run_no_motif_umap.sh      Output: <release>/motif_analysis/no_motif_umap.png, no_motif_umap_clusters.csv
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
MOTIF_DIR="${MOTIF_DIR:-auto}"
[ "$MOTIF_DIR" = auto ] && MOTIF_DIR="$(python ../workflow/release_dir.py)/motif_analysis"
python no_motif_umap.py --dir "$MOTIF_DIR"
