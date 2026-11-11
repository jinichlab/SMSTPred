#!/bin/bash
#SBATCH --job-name=streme_var_scan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/streme_var_scan_%j.log
# After all STREME runs: FIMO cross-scan of every variant's 32 motifs (the repo's run_fimo_cross.sh) and the summary of the variability.
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1   # run from this folder
eval "$(conda shell.bash hook)"; conda activate smstpred_fresh || exit 1
while read -r name src seed; do
    echo "=== FIMO cross $name ==="
    MOTIF_DIR="$PWD/runs/$name" bash ../run_fimo_cross.sh || echo "FIMO failed for $name"
done < variants.txt
python summarise.py
