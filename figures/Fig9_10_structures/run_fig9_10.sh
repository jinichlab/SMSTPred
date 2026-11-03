#!/bin/bash
#SBATCH --job-name=fig9_10
#SBATCH --ntasks=1
#SBATCH --output=run_fig9_10_output.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Figs 9 and 10 of a release: choose the structures by rule (select_structures.py), render them headlessly with ChimeraX, lay them out.
#   [RELEASE_DIR=<release>] sbatch figures/Fig9_10_structures/run_fig9_10.sh     (needs <release>/motif_analysis/structure_search/motif_placements.csv: motif_structure_pick.py)
# ChimeraX: $CHIMERAX, else ~/bin/chimerax. Output: <release>/fig9_10/ (selection.csv, panels/, Fig9_Known.png, Fig10_Distal.png).
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
source "$here/../../workflow/stage_env.sh"
set -e
cd "$here"
mkdir -p "$REL/fig9_10"
python - "$REL" <<'PY'
import sys
import pandas as pd
from select_structures import select
rel = sys.argv[1]
pl = pd.read_csv(f"{rel}/motif_analysis/structure_search/motif_placements.csv")
sel = select(pl)
sel.to_csv(f"{rel}/fig9_10/selection.csv", index=False)
print(sel[["figure", "motif", "why", "pdb", "placement", "d_SAM", "d_substrate"]].to_string(index=False))
PY
export CHIMERAX="$(cfg chimerax)"                 # refresh_config.json "tools" (or $CHIMERAX)
[ -n "$CHIMERAX" ] || unset CHIMERAX
python render_panels.py --release "$REL"
python make_panels.py --release "$REL"
