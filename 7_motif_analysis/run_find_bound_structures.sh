#!/bin/bash
#SBATCH --job-name=find_bound_struct
#SBATCH --ntasks=1
#SBATCH --output=find_bound_structures_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00

# Due-diligence for the apo/AlphaFold-only motifs: for each, expand the FIMO-matching
# reps to their full 50%-id clusters and ask whether any cluster MEMBER carrying the
# motif has a substrate/cofactor-bound PDB we could use instead of the apo/AF model.
#
# NOTE: this job makes outbound HTTPS calls (rest.uniprot.org id-mapping + data.rcsb.org).
# It must run on a node/partition with internet access; if compute nodes are firewalled,
# run it on the login node instead:  bash run_find_bound_structures.sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

python find_bound_structures.py
