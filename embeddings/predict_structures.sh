#!/bin/bash
#SBATCH --job-name=predict_structures
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=predict_structures_output.log

# Predict, with ColabFold, the structures that AlphaFold DB cannot supply for the structure-based embedders (ESM-IF, SaProt): the proteins in
# <structure store>/to_predict.fasta, written by fetch_structures.py (no AlphaFold DB model, or a model for an older UniProt sequence). Then the best model
# of each protein is checked against the current UniProt sequence and stored as <Entry>.pdb (fetch_structures.py --ingest-predictions).
#   python fetch_structures.py --reuse <archived sam_structures>     # first: build the store, write to_predict.fasta
#   sbatch predict_structures.sh                                      # then this (about 3 minutes per protein on one GPU)
#   COLABFOLD_BATCH=/path/to/localcolabfold/colabfold-conda/bin/colabfold_batch sbatch predict_structures.sh
#
# The settings are those of the 26 structures of the archived SAM set (data/README.md): ColabFold 1.6.1, AlphaFold2-ptm, five models, three recycles, no templates,
# no relaxation, the highest-pLDDT model kept. They are passed explicitly, so a change of ColabFold's defaults cannot change them.
# REQUIRES ColabFold, which this repository does not install: the localcolabfold distribution (ColabFold 1.6.1, JAX for CUDA 12, an NVIDIA GPU). How to install it, step by step:
# embeddings/README.md, "Installing ColabFold". Point COLABFOLD_BATCH at its colabfold_batch (or put it on the PATH); the script stops with a message if it is not found.
# NETWORK AND PRIVACY: the MSAs come from the public ColabFold MMseqs2 server (api.colabfold.com), so the sequences of the proteins in to_predict.fasta leave this machine.
# They are public UniProt sequences. No local sequence database is needed.
set -e

here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi

pull="${UNIPROT_PULL_DIR:-$here/../data/uniprot_pull}"
store="${STRUCT_DIR:-$(ls -d "$pull"/release_* | grep -v '/_[^/]*$' | sort | tail -n 1)/structures}"
fasta="$store/to_predict.fasta"
[ -f "$fasta" ] || { echo "$fasta not found: run fetch_structures.py first" >&2; exit 1; }
n=$(grep -c '>' "$fasta" || true)
if [ "$n" -eq 0 ]; then echo "nothing to predict"; exit 0; fi
echo "predicting $n structures from $fasta"

colabfold="$(python "$here/../workflow/cfgget.py" colabfold_batch 2> /dev/null)"; colabfold="${colabfold:-colabfold_batch}"      # $COLABFOLD_BATCH, refresh_config.local.json, refresh_config.json, <repo>/tools/localcolabfold
command -v "$colabfold" > /dev/null || [ -x "$colabfold" ] || { echo "colabfold_batch not found: set COLABFOLD_BATCH to its path (localcolabfold)" >&2; exit 1; }
"$colabfold" --model-type alphafold2_ptm --num-models 5 --num-recycle 3 --rank plddt --msa-mode mmseqs2_uniref_env "$fasta" "$store/colabfold_output"
python -u fetch_structures.py --ingest-predictions "$store/colabfold_output"
# what was run, for the record (ColabFold's server databases change over time, so a prediction is not bit-reproducible: keep the PDBs, they are the record)
{ echo "date $(date -Is)"
  echo "settings: alphafold2_ptm, 5 models, 3 recycles, rank plddt, msa-mode mmseqs2_uniref_env (public ColabFold server), no templates, no relaxation"
  echo "colabfold $("$(dirname "$(command -v "$colabfold")")"/python -m pip show colabfold 2> /dev/null | sed -n 's/^Version: //p')"
  echo "proteins $n"; } > "$store/colabfold_run.txt" || true
