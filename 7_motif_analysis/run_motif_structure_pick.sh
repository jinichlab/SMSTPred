#!/bin/bash
#SBATCH --job-name=motif_struct_pick
#SBATCH --ntasks=1
#SBATCH --output=motif_structure_pick_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
# Step 2 of the motif-structure system: pick the representing structure of every motif by rule and measure it (downloads the mmCIF files it tries from files.rcsb.org; network).
# Needs structure_search/motif_structure_candidates.csv (run_motif_structure_search.sh). Then builds Table 2 (motif_table.csv/.png). Output: <release>/motif_analysis/structure_search/ (motif_placements.csv, motif_pick_log.csv, structures/).
#   [MOTIF_DIR=<release>/motif_analysis | auto] sbatch run_motif_structure_pick.sh
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
python motif_structure_pick.py --dir "$MOTIF_DIR"
python build_motif_table.py --dir "$MOTIF_DIR"      # Table 2: motifs of at least 2-fold, with fold enrichment, placement and distances (motif_table.csv/.png)
