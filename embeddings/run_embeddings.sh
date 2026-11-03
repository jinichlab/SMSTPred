#!/bin/bash
#SBATCH --job-name=embed_new
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=run_embeddings_output.log

# Add the new and changed proteins of a pull to the embedding files. Incremental: a file that already holds a protein's vector for the same sequence (or the same structure
# file) is left alone; only what is new or changed is embedded (incremental.py), so a run on an up-to-date file costs seconds.
#   sbatch run_embeddings.sh                       # sam (default): ProtT5 + ProstT5 of the SAM set (in_deployment, <= 2,000 aa)
#   EMBED_WHAT=up sbatch run_embeddings.sh         # ProtT5 + ProstT5 of the UP set
#   EMBED_WHAT=up-benchmark sbatch run_embeddings.sh   # ProtT5 only of the UP BENCHMARK set (<= 1,000 aa, the proteins mmseqs_clusters.py up clustered): what the SAM-vs-UP comparison trains on
#   EMBED_WHAT=proteomes sbatch run_embeddings.sh      # ProtT5 of the reference proteomes' methyltransferase-domain proteins (IPR029063, <= 2,000 aa; ~116,000), then one H5 per proteome
#   EMBED_WHAT=structures sbatch run_embeddings.sh # SaProt + ESM-IF of the SAM benchmark set, from the structure store of fetch_structures.py
#   EMBED_WHAT=all sbatch run_embeddings.sh
# Where things are (environment variables, defaults in brackets):
#   UNIPROT_PULL_DIR [<repo>/data/uniprot_pull]  the pull; the working embedding files go to <release>/embeddings/
#   ARCHIVED_EMBEDDINGS [$SMSTPRED_DATA or <repo>/data]  the archived embeddings_*.h5 (Zenodo). A working file is first copied from there (cp -n: never over an existing one),
#       then extended; the archived file is only read.
#   ARCHIVED_SEQUENCES [$ARCHIVED_EMBEDDINGS/labels_SAM_full.csv ; UP: $UP_OLD_FASTA, /home/nroethler/Code/files/UP_swiss.fasta]  what the archived vectors were made from (the ledger)
#   ARCHIVED_STRUCTURES  folder of the archived <Entry>.pdb (the ledger of the structure embedders)
#   EMBED_ENV [smstpred_embed] (ProtT5, ProstT5, SaProt), ESMIF_ENV [smstpred_esmif] (ESM-IF): the conda environments of `make embeddings-env`
set -e

here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
[ -z "${FOLDSEEK_BIN:-}" ] && FOLDSEEK_BIN="$(python3 -c "import json;print(json.load(open('$here/../refresh_config.json')).get('tools',{}).get('foldseek_bin',''))" 2> /dev/null)" && [ -n "$FOLDSEEK_BIN" ] && export FOLDSEEK_BIN      # refresh_config.json "tools"
command -v conda > /dev/null && eval "$(conda shell.bash hook)"
use() { command -v conda > /dev/null && conda activate "$1" || { echo "conda environment $1 not found: run 'make embeddings-env' first" >&2; exit 1; }; }

pull="${UNIPROT_PULL_DIR:-$here/../data/uniprot_pull}"
rel="$(ls -d "$pull"/release_* | grep -v '/_' | sort | tail -n 1)"
arch="${ARCHIVED_EMBEDDINGS:-${SMSTPRED_DATA:-$here/../data}}"
work="$rel/embeddings"; mkdir -p "$work"
inputs="$rel/embedding_inputs"
what="${EMBED_WHAT:-sam}"
echo "release $rel | archived files $arch | working files $work | $what"

# $1 = working file name, $2.. = extra arguments for the embedder; copies the archived file first when there is one
start_file() { [ -f "$work/$1" ] || { [ -f "$arch/$1" ] && cp -n "$arch/$1" "$work/$1" && echo "started $1 from the archived file"; } || true; }

sequences() {  # $1 = sam | up
    use "${EMBED_ENV:-smstpred_embed}"
    python -u make_embedding_inputs.py
    if [ "$1" = sam ]; then ledger="${ARCHIVED_SEQUENCES:-$arch/labels_SAM_full.csv}"; else ledger="${ARCHIVED_SEQUENCES:-${UP_OLD_FASTA:-/home/nroethler/Code/files/UP_swiss.fasta}}"; fi
    S="$(echo "$1" | tr a-z A-Z)"
    fasta="$inputs/$1_deployment.fasta"
    start_file "embeddings_${S}_ProtT5.h5"
    python -u embed_prott5.py  --fasta "$fasta" --out "$work/embeddings_${S}_ProtT5.h5"  --ledger "$ledger"
    start_file "embeddings_${S}_ProstT5.h5"
    python -u embed_prostt5.py --fasta "$fasta" --out "$work/embeddings_${S}_ProstT5.h5" --no-per-residue --ledger "$ledger"
}

up_benchmark() {   # ProtT5 of the UP benchmark proteins into embeddings_UP_ProtT5.h5 (started from the archived file; only new / changed proteins are embedded)
    use "${EMBED_ENV:-smstpred_embed}"
    python -u make_embedding_inputs.py --flag in_benchmark
    ledger="${ARCHIVED_SEQUENCES:-${UP_OLD_FASTA:-/home/nroethler/Code/files/UP_swiss.fasta}}"
    start_file embeddings_UP_ProtT5.h5
    python -u embed_prott5.py --fasta "$inputs/up_benchmark.fasta" --out "$work/embeddings_UP_ProtT5.h5" --ledger "$ledger"
}

proteomes() {   # ProtT5 (the training pipeline) of the IPR029063 proteins of the reference proteomes, then split per proteome for 3_predictions/predict.py
    use "${EMBED_ENV:-smstpred_embed}"
    python -u make_proteome_inputs.py --out "$inputs"
    python -u embed_prott5.py --fasta "$inputs/proteomes_ipr029063.fasta" --out "$work/embeddings_proteomes_ProtT5.h5"
    python -u split_proteome_embeddings.py
}

structures() {
    store="$rel/structures"
    [ -d "$store" ] || { echo "$store not found: run fetch_structures.py first" >&2; exit 1; }
    ledger="${ARCHIVED_STRUCTURES:-$arch/sam_structures}"
    [ -d "$ledger" ] || { echo "ARCHIVED_STRUCTURES ($ledger) not found: extract sam_structures.tar.gz there or set it" >&2; exit 1; }
    use "${EMBED_ENV:-smstpred_embed}"
    start_file embeddings_SAM_SaProt.h5
    python -u embed_saprot.py --struct-dir "$store" --out "$work/embeddings_SAM_SaProt.h5" --ledger "$ledger"
    use "${ESMIF_ENV:-smstpred_esmif}"
    start_file embeddings_SAM_ESM-IF.h5
    python -u embed_esm_if_local.py --struct-dir "$store" --pattern '*.pdb' --out "$work/embeddings_SAM_ESM-IF.h5" --ledger "$ledger"
}

case "$what" in
    sam) sequences sam ;;
    up) sequences up ;;
    up-benchmark) up_benchmark ;;
    proteomes) proteomes ;;
    structures) structures ;;
    all) sequences sam; sequences up; structures ;;
    *) echo "EMBED_WHAT must be sam | up | up-benchmark | proteomes | structures | all" >&2; exit 1 ;;
esac
