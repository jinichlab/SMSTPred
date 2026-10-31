#!/bin/bash
#SBATCH --job-name=fimo_all
#SBATCH --ntasks=1
#SBATCH --output=fimo_all_output.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00

# Real per-protein motif status: scan EVERY class member (not just the 50%-id
# representatives, and not homology-propagated) against its class's 8 STREME
# motifs with FIMO, same command/threshold as the rep-level run:
#     fimo --oc fimo_DNA --thresh 1e-4 streme_DNA/streme.txt reps_DNA.clean.fasta
# Only the target FASTA changes: reps -> all members (class_{cls}.fasta).
# FIMO is local (no network) so this runs fine on a compute node.
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

# 1) sanitize members to the 20 standard AAs (FIMO/STREME reject U/O/B/Z/X-source)
for cls in DNA RNA Protein Metabolite; do
    sed '/^>/! s/[^ACDEFGHIKLMNPQRSTVWY]/X/g' "class_${cls}.fasta" > "class_${cls}.clean.fasta"
done

# 2) FIMO per class: that class's 8 motifs vs ALL its members, p<1e-4
for cls in DNA RNA Protein Metabolite; do
    echo "=== FIMO $cls: $(grep -c '>' class_${cls}.clean.fasta) members vs 8 motifs ==="
    fimo --oc "fimo_all_${cls}" --thresh 1e-4 --verbosity 1 \
        "streme_${cls}/streme.txt" "class_${cls}.clean.fasta"
done

# 3) assemble the real per-protein motif-status dataset
python fimo_all_manifest.py

echo "fimo_all done"
