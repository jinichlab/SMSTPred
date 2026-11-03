# Tools that `make env` does not install

`make env` builds the `smstpred` conda environment (Python, PyTorch, scikit-learn, UMAP, RDKit, and from bioconda the MEME suite, MMseqs2, HMMER and MAFFT).
The refresh pipeline (`workflow/run_refresh.sh`) also calls the six tools below. **`bash workflow/setup_tools.sh` checks them all** (it changes nothing) and says what is
missing; `bash workflow/setup_tools.sh --install foldseek|clean|horizyn` installs those three. Where the pipeline looks for each tool is `refresh_config.json`
(`"tools"` and `"annotation"`); an environment variable of the same name in capitals overrides it for one run.

| Tool | Used by | Why it is not in `environment.yml` | Where the pipeline looks |
|---|---|---|---|
| Foldseek | the bake-off structure baseline; SaProt embeddings | a static binary, not a conda package of ours | `$FOLDSEEK_BIN`, `tools.foldseek_bin`, `<repo>/tools/foldseek/bin`, the PATH |
| ColabFold (localcolabfold) | `embeddings/predict_structures.sh`: the proteins AlphaFold DB cannot supply | its own JAX / CUDA conda environment and AlphaFold2 weights | `$COLABFOLD_BATCH`, `tools.colabfold_batch` |
| CLEAN | Stage 6: EC number of every dark protein | needs its own old PyTorch / ESM versions | `annotation.clean_app_dir`, `annotation.clean_env` |
| Horizyn | Stage 6: reactions scored for every dark protein | its own environment and checkpoints | `annotation.horizyn_dir`, `annotation.horizyn_env` |
| ChimeraX | Figs 9 and 10: the structure panels | a vendor application | `$CHIMERAX`, `tools.chimerax` |
| (SLURM, GPU nodes, outbound internet) | the whole chain | the cluster, not software | see "Network" below |

## Foldseek
`make tools` downloads the static binary (`foldseek-linux-avx2.tar.gz` from mmseqs.com, a few tens of MB) into `tools/foldseek/bin/`. Or install it from bioconda into any environment and put it on the PATH.
Check: `foldseek version`.

## ColabFold
Only needed when `<release>/structures/to_predict.fasta` is not empty (no AlphaFold DB model, or one built on an old sequence). The chain runs `embeddings/predict_structures.sh` right after the
structure fetch; it exits in a second when there is nothing to predict (release 2026_03: nothing was needed). Install the [localcolabfold](https://github.com/YoshitakaMo/localcolabfold)
distribution (ColabFold 1.6.1, an NVIDIA GPU), step by step in `embeddings/README.md`, "Installing ColabFold"; then set `tools.colabfold_batch` to its `colabfold-conda/bin/colabfold_batch`.
The sequences are sent to the public ColabFold MSA server (api.colabfold.com). Every verified structure, AlphaFold DB or ColabFold, is kept in `<pull folder>/_structure_cache/` and reused by
later releases where the sequence is unchanged, so a protein is predicted once.

## CLEAN (Yu et al., Science 2023)
1. `git clone https://github.com/tttianhao/CLEAN.git tools/CLEAN` (the pretrained weights for the 70 % and 100 % identity splits are in the repository, `app/data/pretrained`). Its licence is a non-exclusive research-use licence: read it.
2. `conda env create -f 6_dark_cluster_annotation/clean/environment_clean.yml` (the versions the runs here used).
3. CLEAN downloads ESM-1b (about 2.5 GB) on its first run; the vendored `app/esm` folder has to be importable (`run_clean_dark.sh` puts it on `PYTHONPATH`).
4. Set `annotation.clean_app_dir` to the `app/` folder.
Files the pipeline expects in `app/`: `CLEAN_infer_fasta.py`, `data/inputs/`, `data/pretrained/`; it copies `6_dark_cluster_annotation/clean/CLEAN_infer_fasta_pvalue.py` there itself.
Note: this was written from the checkout used on this machine; a fresh clone has not been tested here.

## Horizyn (Dayhoff Labs)
1. `git clone https://github.com/dayhofflabs/horizyn.git tools/horizyn`.
2. `conda env create -f 6_dark_cluster_annotation/horizyn/environment_horizyn.yml`.
3. The two official checkpoints (`horizyn_v1_0_dev.ckpt`, `horizyn_v1_0_inf.ckpt`, about 201 MB each) from Zenodo record 20348783 into `checkpoints/` (md5 `5b1f938f…` and `cf6775b7…`); `--install horizyn` does this.
4. Set `annotation.horizyn_dir` to the folder holding `checkpoints/`, `configs/` and `horizyn/`. The pipeline never writes into it: `run_horizyn_release.sh` links it into `<release>/dark_annotation/horizyn_run/`.
Note: untested on a fresh clone; the package layout of this machine's checkout is what the scripts assume.

## ChimeraX
Install UCSF ChimeraX (1.11.1 was used) from https://www.cgl.ucsf.edu/chimerax/download.html; a rootless unpack is enough. Headless rendering needs `libGLU.so.1` and `libOpenGL.so.0`; if the node lacks them, a
small launcher script that adds a folder with those two libraries to `LD_LIBRARY_PATH` does it (the one here: `~/bin/chimerax`). Set `tools.chimerax` to the launcher. Check:
`chimerax --offscreen --nogui --exit --cmd "version"`.

## Network
Stages that need outbound HTTPS (from the compute node that runs them): the UniProt pull and the notes (rest.uniprot.org, sparql.rhea-db.org, ftp.expasy.org), the structure fetch
(alphafold.ebi.ac.uk), ColabFold (api.colabfold.com), the structure search and download (data.rcsb.org, files.rcsb.org). If compute nodes are firewalled, run those stages on a node that has access:
the chain scripts are plain bash, `bash <script>` works on a login node.
