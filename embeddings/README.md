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

## Structures for the refresh (ESM-IF, SaProt): one store, built by code

The structure-based embedders run on the benchmark set only (<= 1,000 aa). Both read the same folder of `<Entry>.pdb` files, so they cannot see different structures and the exact
files that were embedded go to Zenodo (`sam_structures.tar.gz`, about 350 MB). `embed_esm_if.py` (streams AlphaFold DB, embeds, deletes) is how the archived ESM-IF file was made and
stays for that; the refresh embeds ESM-IF from the store with `embed_esm_if_local.py --struct-dir <store> --pattern '*.pdb'`.

```bash
python fetch_structures.py --reuse <folder of the archived sam_structures>     # = make structures REUSE=...
sbatch predict_structures.sh                                                    # ColabFold for the ones AlphaFold DB cannot supply
```

`fetch_structures.py` keeps or reuses a structure only if its sequence equals the **current** UniProt sequence, asks the AlphaFold DB API (`alphafold.ebi.ac.uk/api/prediction`, using the
file URL of each record: no model version is hard-coded) for the rest, and accepts a model only if the model's own sequence matches. A protein with no model, or a model built on an older
sequence (UniProt changed it afterwards; an unchecked fetch would embed the structure of a different protein), goes to `to_predict.fasta`. `structures_status.csv` records for every protein
`status` (`kept`, `reused_archived`, `afdb`, `colabfold`, `no_afdb`, `afdb_sequence_mismatch`, `error`), the source, the model version, the pLDDT and a note.
`predict_structures.sh` runs ColabFold with the settings of the archived 26 predictions (AlphaFold2-ptm, five models, three recycles, no templates, no relaxation, best pLDDT; **the
sequences go to the public ColabFold MSA server**) and ingests the rank_001 models through `fetch_structures.py --ingest-predictions`, again only where the sequence matches.
Release 2026_03: 7,553 benchmark proteins, 7,399 reused from the archive, 143 from AlphaFold DB, 11 to predict (7 without a model, 4 with a model for another sequence).

### Installing ColabFold (only needed to predict structures; not installed by `make env` or `make embeddings-env`)

`predict_structures.sh` calls `colabfold_batch`, which this repository does not install. It is the [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) distribution of ColabFold, and you need
all of the following before the script can run:

* a Linux machine with an **NVIDIA GPU and a CUDA 12 driver** (the installer pulls `jax[cuda12]==0.5.3`), `wget` and `git`, and several GB of free disk (its own conda environment plus the AlphaFold2 weights);
* the installer of that project, run from a folder where you want ColabFold to live. It creates `./localcolabfold/` with a private conda environment (Python 3.10, OpenMM, kalign2, HH-suite, MMseqs2),
  installs **ColabFold 1.6.1** and JAX, and downloads the AlphaFold2 weights (`python -m colabfold.download`):
  ```bash
  wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh      # or copy it from the localcolabfold repository
  bash install_colabbatch_linux.sh                                                                          # takes a while; its README describes the current steps
  export COLABFOLD_BATCH="$PWD/localcolabfold/colabfold-conda/bin/colabfold_batch"                          # what predict_structures.sh looks for (or put that bin/ on your PATH)
  ```
* network access to the public ColabFold MSA server at run time (the sequences are sent there, see above); no local sequence database is needed.

The archived 26 structures were made with ColabFold 1.6.1 from this installer; `predict_structures.sh` writes `colabfold_run.txt` (date, settings, ColabFold version) next to the predictions.
Without ColabFold you can still use any other predictor: put its models into a folder and run `fetch_structures.py --ingest-predictions <folder>`; a structure is accepted only if its sequence matches.

**What works for someone else.** Re-running the manuscript needs no structure prediction: the PDBs of the archived set, with the 26 predicted ones, are in `sam_structures.tar.gz` on Zenodo, and `fetch_structures.py --reuse`
treats them as already done. Only a refresh onto a newer UniProt release produces proteins that need a prediction, and `predict_structures.sh` then needs what the repository does not install: ColabFold
(`COLABFOLD_BATCH`, a GPU, a few GB of weights) and network access to the public ColabFold MSA server, a shared service with fair-use limits whose sequence databases change over time, so a prediction made later is not
bit-identical to an earlier one. That is why the PDBs are archived instead of regenerated, and why `fetch_structures.py --ingest-predictions` accepts a folder from any run: it keeps a structure only if its sequence matches.


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

## `EMBED_WHAT=up-benchmark` (the SAM-vs-UP comparison)

`sbatch run_embeddings.sh` with `EMBED_WHAT=up-benchmark` embeds **ProtT5 only** for the UP benchmark set (`in_benchmark`, <= 1,000 aa, the proteins `mmseqs_clusters.py up` clustered) into
`<release>/embeddings/embeddings_UP_ProtT5.h5`, started from the archived UP file, so only new or changed proteins are embedded. Release 2026_03: 204,941 proteins, of which 200,509 reused and 4,432 (4,302 new,
130 with a changed sequence) embedded in 3.5 minutes. `EMBED_WHAT=up` (ProtT5 + ProstT5 of the deployment set) is unchanged.

## `EMBED_WHAT=proteomes` (the reference proteomes, ProtT5)

`make_proteome_inputs.py` selects the proteins with the methyltransferase domain InterPro IPR029063 and <= 2,000 aa from `<release>/reference_proteomes/*.csv` (the longer ones, from `reference_proteomes_long/`, are
listed in `proteomes_not_predicted.csv` and get no prediction: the longest training protein is 1,848 aa); `embed_prott5.py` embeds them with the same pipeline as the training vectors into
`embeddings_proteomes_ProtT5.h5`; `split_proteome_embeddings.py` writes one `<proteome>_embeddings.h5` per proteome (`<release>/embeddings/proteomes_ProtT5/`), the layout of the archived tree-of-life files, which
`3_predictions/predict.py` scores one file at a time. Release 2026_03: 115,850 proteins in 1,249 proteomes (1,703 longer than 2,000 aa not predicted), 62 minutes on one GPU. The archived tree-of-life vectors were not
reused (no sequences to check them against); they agree with the training pipeline anyway (2_training/README.md).
