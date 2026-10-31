# Embeddings: how the archived vectors were made

Everything downstream (clustering, training, prediction) starts from per-protein embeddings. The embeddings themselves are archived on
Zenodo (see [`../data/README.md`](../data/README.md)); the scripts here re-compute them. Re-computing is optional: the analyses read the archived files.

| Archived file(s) | Embedder (this folder) | Model | Vector | Details |
|---|---|---|---|---|
| `embeddings_SAM_ProtT5.h5`, `embeddings_UP_ProtT5.h5` | `embed_prott5.py` (defaults) | `Rostlab/prot_t5_xl_half_uniref50-enc` | 1024, float16 | mean over the L residue tokens, fp16 model on a GPU |
| `embeddings_tree_of_life_ProtT5.tar` | `embed_prott5.py --model Rostlab/prot_t5_xl_uniref50 --pool masked-mean --precision full` | `Rostlab/prot_t5_xl_uniref50` | 1024, float32 | mean over every non-padding token (end-of-sequence token included), fp32 |
| `embeddings_SAM_ProstT5.h5`, `embeddings_UP_ProstT5.h5`, `embeddings_tree_of_life_ProstT5.tar` | `embed_prostt5.py --no-per-residue` | `Rostlab/ProstT5` | 1024, float16 | `<AA2fold>` prefix, encoder only, mean over the residue tokens, fp16 |
| `embeddings_SAM_ESM-IF.h5` | `embed_esm_if.py` (AlphaFold DB models, version v6) and `embed_esm_if_local.py` (the 26 ColabFold structures) | ESM-IF1 (`esm_if1_gvp4_t16_142M_UR50`), fair-esm 2.0.0 | 512, float32 | chain A, mean over residues of the encoder output |
| `embeddings_SAM_SaProt.h5` | `embed_saprot.py` | `westlake-repl/SaProt_650M_AF2` | 1280, float32 | structure-aware sequence (amino acid + Foldseek 3Di letter per residue), mean over residues |

Input sequences were limited to 1,000 residues or fewer (filter your FASTA first); non-standard residues (U, Z, O, B) are mapped to X.
The ProstT5 files of the UP and tree-of-life sets were written by batched variants of `embed_prostt5.py` that produce the same vector.
ProtT5 of the training and UP sets was computed with the half-precision encoder of the ProtT5 authors' reference embedder; `embed_prott5.py` is an
independent implementation of the same procedure (its output is checked against the archived vectors below).

## Environments

Two, because fair-esm 2.0.0 needs the PyTorch-Geometric extensions built for one torch release:

```bash
make embeddings-env                  # smstpred_embed (ProtT5, ProstT5, SaProt) and smstpred_esmif (ESM-IF), from the two environment_*.yml here
make tools                           # Foldseek binary for SaProt (or point FOLDSEEK_BIN at your own)
```

Versions used for the archived files: transformers 4.57.1, sentencepiece 0.2.1, h5py 3.15, torch nightly 2.12 (CUDA 12.8) for the ProtT5 / ProstT5 /
SaProt runs; fair-esm 2.0.0, torch 2.6.0 (CUDA 12.4), biotite 0.39, torch-geometric 2.7.0 for ESM-IF; Foldseek build 718d4217. The environment
files pin released versions (torch 2.8.0 instead of the nightly).

## Re-computing

```bash
python embed_prott5.py  --fasta proteins.fasta --out prott5.h5
python embed_prostt5.py --fasta proteins.fasta --out prostt5.h5 --no-per-residue
python embed_saprot.py  --struct-dir ../data/sam_structures --out saprot.h5      # <Entry>.pdb files; needs Foldseek
python embed_esm_if.py  --csv ../data/labels_SAM_full.csv --out esmif.h5         # downloads AlphaFold DB models (needs network)
python embed_esm_if_local.py --struct-dir ../data/sam_structures --pattern '*.pdb' --out esmif_local.h5
```

Output layouts: ProtT5, ESM-IF and SaProt write one dataset per accession; ProstT5 writes `/<accession>/mean` (and `/<accession>/residue` unless
`--no-per-residue`). The sequences of the SAM set are in `../data/labels_SAM_full.csv`, the structures in `../data/sam_structures/` (AlphaFold DB models
plus the 26 locally predicted ones).

## Check against the archive

`make embeddings-check` (or `sbatch embeddings/run_spotcheck.sh`; one GPU, about 30 minutes once the environments exist, mostly loading the models; the first run also downloads
them from Hugging Face, several GB) re-embeds 20 SAM proteins (seed 0) with every embedder and compares each new vector with the archived one.

Result of the last run (RTX 6000 Ada, the two environments built from the `environment_*.yml` files here, Foldseek build 718d4217):

| Embedder | Minimum cosine similarity to the archived vector (20 proteins) | Largest absolute difference in any dimension |
|---|---|---|
| ProtT5 (`embed_prott5.py`) | 1.000000 | 1.2e-4 |
| ProstT5 (`embed_prostt5.py`) | 1.000000 | 4.9e-4 |
| SaProt (`embed_saprot.py`) | 1.000000 | 1.8e-7 |
| ESM-IF (`embed_esm_if_local.py`, on the archived structures) | 1.000000 | 2.7e-7 |

The float16 files (ProtT5, ProstT5) differ at half-precision rounding, the float32 files at single-precision rounding. The script exits with an error if any minimum cosine falls
below 0.99 (`--min-cosine`). The report is written to `spotcheck_output/` (`spotcheck_report_<models>.csv`, one row per model and protein).
