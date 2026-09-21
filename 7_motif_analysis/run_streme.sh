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
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

# 1) sanitize reps to the 20 standard AAs (STREME rejects U/O/etc)
for cls in Other Protein RNA DNA; do
    sed '/^>/! s/[^ACDEFGHIKLMNPQRSTVWY]/X/g' "reps_${cls}.fasta" > "reps_${cls}.clean.fasta"
done

# 2) build each class's control = concatenation of the other three classes' reps
for cls in Other Protein RNA DNA; do
    : > "ctrl_${cls}.fasta"
    for o in Other Protein RNA DNA; do
        [ "$o" = "$cls" ] && continue
        cat "reps_${o}.clean.fasta" >> "ctrl_${cls}.fasta"
    done
done

# 3) STREME per class -- full positives, full control, width capped 8-18
for cls in Other Protein RNA DNA; do
    echo "=== STREME $cls: $(grep -c '>' reps_${cls}.clean.fasta) pos vs $(grep -c '>' ctrl_${cls}.fasta) ctrl ==="
    streme --p "reps_${cls}.clean.fasta" --n "ctrl_${cls}.fasta" --protein \
        --minw 8 --maxw 18 --nmotifs 8 --oc "streme_${cls}"
done
echo "STREME done -- logos in streme_{Other,Protein,RNA,DNA}/streme.html"
