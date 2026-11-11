#!/bin/bash
#SBATCH --job-name=distal_contacts
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=run_%j.log
# REL_NEW: release folder of the rebuilt data; REL_OLD: the archived release folder (both hold motif_analysis/structure_search/)
: "${REL_NEW:?set REL_NEW}"; : "${REL_OLD:?set REL_OLD}"
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source ~/.bashrc
python distal_motif_contacts.py --dir "$REL_NEW/motif_analysis" --motifs DNA-4 Metabolite-7 Protein-6 RNA-5 RNA-7 --out out_rebuilt --cif-cache cif
python distal_motif_contacts.py --dir "$REL_OLD/motif_analysis" --motifs DNA-4 Metabolite-3 Protein-6 RNA-6 RNA-8 --out out_archived --cif-cache cif
