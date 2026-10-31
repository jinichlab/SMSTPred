#!/bin/bash
#SBATCH --job-name=embed_spotcheck
#SBATCH --ntasks=1
#SBATCH --output=embeddings_spotcheck_output.log
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00

# Re-embed 20 SAM proteins with every embedder of this folder and compare with the archived vectors (needs the Zenodo data: `make data`).
# Two environments are used, created from the two environment_*.yml files here (`make embeddings-env`): ProtT5, ProstT5 and SaProt run in
# smstpred_embed, ESM-IF in smstpred_esmif (rename them with EMBED_ENV / ESMIF_ENV). SaProt needs Foldseek (FOLDSEEK_BIN, PATH, or `make tools`).
# Submit as `sbatch embeddings/run_spotcheck.sh` from anywhere, or run `bash embeddings/run_spotcheck.sh` on a GPU machine.
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1

use_env() {      # use_env <name>: activate it unless it is the active environment already
    if [ "${CONDA_DEFAULT_ENV:-}" != "$1" ] && command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "$1" || { echo "conda environment $1 not found: create it with 'make embeddings-env'" >&2; exit 1; }
    fi
}

status=0
use_env "${EMBED_ENV:-smstpred_embed}"
python spotcheck_embeddings.py --models prott5 prostt5 saprot || status=1
use_env "${ESMIF_ENV:-smstpred_esmif}"
python spotcheck_embeddings.py --models esmif || status=1
exit $status
