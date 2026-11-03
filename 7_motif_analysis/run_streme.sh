#!/bin/bash
#SBATCH --job-name=motif_streme
#SBATCH --ntasks=1
#SBATCH --output=streme_output.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00

# Discriminative motif discovery: for each substrate class, STREME finds motifs
# enriched in that class vs the others. We use the FULL set of 50%-identity family
# representatives (no random subsampling) so no family is dropped; the only speed
# concession is capping motif width at 8-18 (the wide-motif refinement was what
# made the unbounded run take ~8.5 h/class). Control = all reps of the other three
# classes. One STREME run per class.
# Environment: the conda environment that is already active, else the one `make env` builds (SMSTPRED_ENV renames it).
# Runs as `sbatch run_*.sh` (the job inherits the submitting shell's environment) or as plain `bash run_*.sh`.
# The scripts use relative paths, so start in this script's folder (under sbatch $0 is a spooled copy: ask SLURM where it came from).
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi

# MOTIF_DIR=<release>/motif_analysis runs on a refresh release (after run_dedup.sh made its reps_*.fasta); unset: this folder.
# STREME_PREP=0 skips steps 1-2 and STREME_CLASSES="Protein" (a subset) runs only those classes, so that the four classes can run as parallel jobs
# (STREME is single-threaded): one job with STREME_PREP=1 STREME_CLASSES="" first, then one per class. Unset: everything, one class after the other.
[ "${MOTIF_DIR:-}" = auto ] && MOTIF_DIR="$(python "$PWD/../workflow/release_dir.py")/motif_analysis"      # auto: the newest release (the refresh chain)
[ -n "${MOTIF_DIR:-}" ] && cd "$MOTIF_DIR"
PREP="${STREME_PREP:-1}"
CLASSES="${STREME_CLASSES-Metabolite Protein RNA DNA}"

# 1) sanitize reps to the 20 standard AAs (STREME rejects U/O/etc)
if [ "$PREP" = 1 ]; then
for cls in Metabolite Protein RNA DNA; do
    sed '/^>/! s/[^ACDEFGHIKLMNPQRSTVWY]/X/g' "reps_${cls}.fasta" > "reps_${cls}.clean.fasta"
done
fi

# 2) build each class's control = concatenation of the other three classes' reps
if [ "$PREP" = 1 ]; then
for cls in Metabolite Protein RNA DNA; do
    : > "ctrl_${cls}.fasta"
    for o in Metabolite Protein RNA DNA; do
        [ "$o" = "$cls" ] && continue
        cat "reps_${o}.clean.fasta" >> "ctrl_${cls}.fasta"
    done
done
fi

# 3) STREME per class -- full positives, full control, width capped 8-18
for cls in $CLASSES; do
    echo "=== STREME $cls: $(grep -c '>' reps_${cls}.clean.fasta) pos vs $(grep -c '>' ctrl_${cls}.fasta) ctrl ==="
    streme --p "reps_${cls}.clean.fasta" --n "ctrl_${cls}.fasta" --protein \
        --minw 8 --maxw 18 --nmotifs 8 --oc "streme_${cls}"
done
echo "STREME done -- logos in streme_{Metabolite,Protein,RNA,DNA}/streme.html"
