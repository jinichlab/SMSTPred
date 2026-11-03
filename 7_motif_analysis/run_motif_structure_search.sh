#!/bin/bash
#SBATCH --job-name=motif_struct_search
#SBATCH --ntasks=1
#SBATCH --output=motif_structure_search_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
# Which FIMO carriers of each motif have an experimental structure with SAM and a substrate? (needs outbound HTTPS to rest.uniprot.org and data.rcsb.org: if the compute
# nodes are firewalled run it on the login node: bash run_motif_structure_search.sh)
#   [MOTIF_DIR=<release>/motif_analysis | auto] sbatch run_motif_structure_search.sh      Output: <release>/motif_analysis/structure_search/
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
python motif_structure_search.py --dir "$MOTIF_DIR"
