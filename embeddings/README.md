# Embeddings

ProtT5-XL-U50 protein language model embedding generation pipeline for SAM-dependent methyltransferases.

## Scripts

| Script | Description |
|--------|-------------|
| `batched_prott5_multiple_files.py` | Batch embedding generation across multiple reference proteomes. Filters for IPR029063 domain (SAM MTase marker), proteins ≤1000 aa, and generates 1024-dimensional mean-pooled embeddings. |
| `prott5_embedder.py` | General-purpose ProtT5 embedder. Converts FASTA sequences to HDF5 embeddings with configurable per-protein or per-residue output. |
| `create_fasta.py` | Converts CSV files (with Entry, Sequence columns) to FASTA format for embedding input. |
| `match_embeddings.py` | Utility to merge HDF5 embeddings with DataFrames by protein Entry ID. |
| `script_prott5.sh` | Shell wrapper for running the embedding pipeline. |

## Usage

```bash
# Generate embeddings from FASTA
python prott5_embedder.py -i sequences.fasta -o embeddings.h5 --per_protein 1

# Batch process multiple proteomes
python batched_prott5_multiple_files.py
```

## Input/Output

- **Input**: CSV files with protein sequences and InterPro annotations, or FASTA files
- **Output**: HDF5 files containing 1024-dimensional per-protein embeddings

## Environment

Uses the module-specific `environment.yml` in this directory, or the root-level environment.
