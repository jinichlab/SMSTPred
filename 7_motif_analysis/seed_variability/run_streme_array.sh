#!/bin/bash
#SBATCH --job-name=streme_var
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=0-27
#SBATCH --output=logs/streme_var_%A_%a.log
# One STREME run (one class, one variant) per array task: task = 4 * variant + class. Same settings as 7_motif_analysis/run_streme.sh; only --seed differs.
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1   # run from this folder
eval "$(conda shell.bash hook)"; conda activate smstpred_fresh || exit 1
classes=(Metabolite Protein RNA DNA)
v=$((SLURM_ARRAY_TASK_ID / 4)); cls=${classes[$((SLURM_ARRAY_TASK_ID % 4))]}
read -r name src seed < <(sed -n "$((v + 1))p" variants.txt)
src=$(eval echo "$src"); [ -d "$src" ] || { echo "set ARCHIVED_MOTIF_DIR and REBUILT_MOTIF_DIR (variants.txt) to the two motif_analysis folders"; exit 1; }
mkdir -p "runs/$name"; cd "runs/$name" || exit 1
for f in corpus.csv "reps_${cls}.clean.fasta" "ctrl_${cls}.fasta"; do ln -sf "$src/$f" "$f"; done
echo "=== $name $cls seed $seed: $(grep -c '>' reps_${cls}.clean.fasta) pos vs $(grep -c '>' ctrl_${cls}.fasta) ctrl ==="
streme --p "reps_${cls}.clean.fasta" --n "ctrl_${cls}.fasta" --protein --minw 8 --maxw 18 --nmotifs 8 --seed "$seed" --oc "streme_${cls}"
