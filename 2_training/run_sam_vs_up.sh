#!/bin/bash
#SBATCH --job-name=sam_vs_up
#SBATCH --ntasks=1
#SBATCH --array=1-16
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=run_sam_vs_up_%A_%a.log

# Stage 2: SAM-trained vs UP-trained ProtT5 MLP on the same held-out SAM test clusters, one split per array task (sam_vs_up.py). Needs the bake-off (its MLP_ProtT5 is the SAM arm), the UP splits
# (1_clustering/up_splits.py) and the UP ProtT5 embeddings (embeddings/run_embeddings.sh, EMBED_WHAT=up-benchmark). About 5-10 minutes per split on one GPU.
#   sbatch run_sam_vs_up.sh                     # all 16 splits
#   sbatch --array=1 --export=ALL,SMOKE=1 run_sam_vs_up.sh      # tiny end-to-end test into <release>/sam_vs_up/_smoke
#   SAM_VS_UP_ARMS=SAM,UP sbatch run_sam_vs_up.sh      # without the SAMclean sensitivity arm
# Then (CPU): sbatch run_sam_vs_up_report.sh. Release folder: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR. Resumable.
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
python -u sam_vs_up.py --iteration "${SLURM_ARRAY_TASK_ID:-${ITERATION:?run with sbatch (SLURM_ARRAY_TASK_ID) or set ITERATION}}" ${SMOKE:+--smoke}
