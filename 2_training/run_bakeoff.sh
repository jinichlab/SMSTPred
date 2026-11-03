#!/bin/bash
#SBATCH --job-name=bakeoff
#SBATCH --ntasks=1
#SBATCH --array=1-16
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=run_bakeoff_%A_%a.log

# Stage 2: the classifier bake-off on the 16 cluster-disjoint train / validation / test splits, one split per array task (see bakeoff.py for the protocol: six embedding MLPs trained with
# BCE for 2,000 epochs with the best-validation-epoch weights kept, MMseqs2 / profile-HMM / Foldseek baselines, a fixed 0.5 threshold). About 10 minutes per split on one GPU.
#   sbatch run_bakeoff.sh                       # all 16 splits (needs 1_clustering/tvt_splits.py, the embeddings and the structures of the release)
#   sbatch --array=3 run_bakeoff.sh             # one split
#   sbatch --array=1 --export=ALL,SMOKE=1 run_bakeoff.sh      # tiny end-to-end test (20 epochs, 900 proteins) into <release>/benchmark/_smoke
# Then (CPU):  sbatch run_bakeoff_report.sh
# Release folder: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR (default <repo>/data/uniprot_pull). Settings: ../refresh_config.json. Resumable: a finished split is skipped.
# Environment: the conda environment that is already active, else the one `make env` builds (SMSTPRED_ENV renames it); Foldseek: `make tools` (or $FOLDSEEK_BIN).
set -e

here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
[ -z "${FOLDSEEK_BIN:-}" ] && FOLDSEEK_BIN="$(python3 -c "import json;print(json.load(open('$here/../refresh_config.json')).get('tools',{}).get('foldseek_bin',''))" 2> /dev/null)" && [ -n "$FOLDSEEK_BIN" ] && export FOLDSEEK_BIN      # refresh_config.json "tools"
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
python -u test_bakeoff.py
python -u bakeoff.py --iteration "${SLURM_ARRAY_TASK_ID:-${ITERATION:?run with sbatch (SLURM_ARRAY_TASK_ID) or set ITERATION}}" ${SMOKE:+--smoke}
