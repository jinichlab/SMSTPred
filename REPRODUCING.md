# Reproducing the analyses

One entry point: the `Makefile`. The archived Zenodo record holds the inputs and the results of one UniProt release (2026_03). `make data` rebuilds that release's folder from it,
`make verify` re-runs every stage that works on a CPU from those files in a scratch copy and compares what each stage writes with the archived file, and `make verify-gpu` re-trains the
deployed ensemble and compares outcomes with stated tolerances. The steps are registered in [`workflow/steps.py`](workflow/steps.py); the tables in [section 3](#3-pass--fail-tables) are the output.
(The registry of the 2024-era paper, with the iTOL tree, 0.6 clustering and ProstT5 deployment, is kept as `workflow/legacy_steps_paper_2024.py` and is not run.)

## 1. Quick start

```bash
git clone https://github.com/jinichlab/SMSTPred.git && cd SMSTPred
make env                          # conda environment "smstpred" from environment.yml (Python 3.12, pinned versions; also MEME suite, MMseqs2, HMMER, MAFFT)
conda activate smstpred
make data DATA_SRC=zenodo:10.5281/zenodo.23655998      # download, checksum, and assemble the release folder data/uniprot_pull/release_2026_03/
#   already have the files?  make data DATA_SRC=/path/to/files LINK=--link
make check                        # SHA-256 of every file in data/ against MANIFEST.sha256
make verify                       # CPU: unit tests, clustering, the ensemble applied to every proteome, tree, Figs 2-8 and 11, Tables S1-S2, motif picks; writes verify/report.json
make verify-gpu                   # GPU: re-trains the deployment ensemble and compares outcomes (tolerances below)
```

`sbatch workflow/slurm/verify_cpu.sh` and `sbatch workflow/slurm/verify_gpu.sh` run the same as single SLURM jobs; `workflow/slurm/test_fresh_clone.sh` builds a clean clone and environment first.

| Target | What it does |
|---|---|
| `make env` | creates the conda environment |
| `make data` | puts the Zenodo files in `data/` (or `$SMSTPRED_DATA`) and assembles the release folder (`workflow/fetch_data.py`, the inverse of `workflow/build_zenodo_files.py`) |
| `make check` | verifies the checksums of `data/` |
| `make tools` / `make annotation-test` / `make bakeoff-test` | Foldseek binary; offline unit tests of the stages |
| `make steps` | lists every step: id, stage, command |
| `make stage1` … `make stage7`, `make figures` | runs one stage in place on the release folder |
| `make verify` / `make verify-gpu` | the checks above; `RELEASE_DIR=<folder>` points them at another release |
| `make refresh` | the whole chain for the newest UniProt release (`bash workflow/run_refresh.sh --help`): needs network, a GPU and days; `make refresh-dry-run` prints the job chain |

`make verify` never touches the archived release: it copies the repository (not `data/`) to `.verify_work/` and the release folder to `.verify_work_release/` (links resolved), runs each stage's own driver
there with `RELEASE_DIR` set to the copy, and compares every file a step wrote or changed with the archived file of the same relative path. One step: `python workflow/run_steps.py --mode verify --only s4_fig5`.

## 2. What "reproduces" means

| File type | Comparison | Tolerance |
|---|---|---|
| tables and statistics (`.csv`, `.csv.gz`, `.tsv`, `.txt`, `.md`, `.json`) | same text outside the numbers; every number within tolerance | relative 1e-5 (1e-3 for the ensemble predictions, which are float32) |
| figures (`.png`) | same pixel size (within 2%); mean absolute grey-level difference | at most 12 of 255: catches a missing panel or a changed layout, not anti-aliasing or font hinting |
| vector graphics, HTML, binary arrays and models | not compared | the numbers behind them are compared through the tables |

**PASS**: at least one regenerated file matched its archived version and none differed. **BUILT**: the step ran and wrote output but nothing archived corresponds to it (the unit tests, the clustering). **FAIL**: a step crashed, timed out, or a file differs beyond tolerance.

## 3. Pass / fail tables

### CPU steps (`make verify`)

<!-- VERIFY-TABLE-START -->
Last verification: 2026-11-04T13:17:14 (Python 3.12.15, numpy 2.1.3, pandas 2.2.3, scipy 1.15.2, torch 2.8.0, matplotlib 3.10.9). 13 of 17 steps reproduce a committed result, 4 build their output without a committed reference to compare, 0 fail.

| Step | Reproduces | Result | Files compared | Time |
|---|---|---|---|---|
| `t_binning` | Substrate typing rules (Rhea compound classes -> Metabolite / Protein / RNA / DNA) | built (no reference) | ran; wrote no files | 1 s |
| `t_training` | Cluster-disjoint splits, bake-off, MTase-vs-UP and ensemble code, MCC threshold search | built (no reference) | ran; wrote no files | 60 s |
| `t_annotation` | Cluster selection, dark-protein inputs, CLEAN summary, Horizyn family table, Fig 8 rule, motif structure search and pick | built (no reference) | ran; wrote no files | 5 s |
| `s1_cluster` | 502 clusters of the benchmark set, 16 splits, cluster columns of the label table | built (no reference) | 2 new-file | 14 s |
| `s2_bakeoff_report` | Fig 2A statistics, Table 1, Table S4, the train / validation loss figures | PASS | 20 identical | 16 s |
| `s2_sam_vs_up_report` | Fig 2B: MTase-trained vs UP-trained, Wilcoxon-Holm and Nadeau-Bengio | PASS | 8 identical | 5 s |
| `s2_fig2` | Fig 2 composite (A + B) | PASS | 2 identical | 5 s |
| `s3_predict` | Scores, calls, ensemble disagreement and applicability-domain flag for every proteome protein | PASS | 1 identical, 1249 new-file | 28 s |
| `s4_tree` | Fig 3: taxonomy tree and per-organism annotation | PASS | 11 identical, 2 skipped | 8 s |
| `s4_fig4` | Fig 4 and its Mann-Whitney / Kruskal-Wallis statistics | PASS | 7 identical | 4 s |
| `s4_fig5` | Fig 5 power-law fits, level differences, clade-bootstrap CIs, lineage statistics | PASS | 9 identical, 1 skipped | 19 s |
| `s5_clusters` | Dark / metazoan-protein / vertebrate-DNA cluster selection and Fig 7 | PASS | 11 identical | 6 s |
| `s6_tables` | Horizyn family table, CLEAN summary, Tables S1 and S2, Fig 7 with the Fig 8 clusters | PASS | 17 identical | 17 s |
| `s6_fig8` | Fig 8: Horizyn reaction profiles of the rule-chosen dark clusters | PASS | 2 identical | 9 s |
| `s7_no_motif_umap` | Fig 11: enzymes with no discriminative motif on the map | PASS | 2 identical | 2 s |
| `s7_motif_pick` | Motif structure pick by rule, Table 2 | PASS | 77 identical, 185 new-file | 68 s |
| `s9_key_numbers` | key_numbers.md (every number quoted in the text) and the figure index | PASS | 1 identical, 23 new-file | 4 s |
<!-- VERIFY-TABLE-END -->

### GPU re-run (`make verify-gpu`)

A retrained network does not reproduce to the last digit across PyTorch / CUDA / GPU generations, so the deployment step compares outcomes: the same proteins, clusters and folds, the fixed thresholds, per-label out-of-fold MCC (within 0.02), out-of-fold scores (max 0.1, mean 0.01) and the applicability-domain cutoff (0.01) of the six-member ensemble. The tolerances are provisional (`workflow/gpu_checks.py`) until more GPU generations have been compared; the manuscript's numbers are those of the archived ensemble.

<!-- VERIFY-GPU-TABLE-START -->
Last verification: 2026-11-04T13:17:14 (Python 3.12.15, numpy 2.1.3, pandas 2.2.3, scipy 1.15.2, torch 2.8.0, matplotlib 3.10.9). 1 of 1 steps reproduce a committed result, 0 build their output without a committed reference to compare, 0 fail.

| Step | Reproduces | Result | Checks | Time |
|---|---|---|---|---|
| `b_deploy` | Deployment ensemble re-trained from the labels: thresholds, out-of-fold MCC, out-of-fold scores | PASS | 8 of 8 passed | 54 s |
| &nbsp;&nbsp;`b_deploy` check | same proteins, clusters and folds | PASS | 7582 proteins | |
| &nbsp;&nbsp;`b_deploy` check | six ensemble members written | PASS |  | |
| &nbsp;&nbsp;`b_deploy` check | every label called at the same fixed threshold as archived | PASS | {'Metabolite': 0.5, 'Protein': 0.5, 'RNA': 0.5, 'DNA': 0.5} | |
| &nbsp;&nbsp;`b_deploy` check | overall out-of-fold MCC within 0.01 | PASS | 0.9551 vs 0.9551 | |
| &nbsp;&nbsp;`b_deploy` check | per-label out-of-fold MCC within 0.02 | PASS | Metabolite 0.962 vs 0.962, Protein 0.913 vs 0.913, RNA 0.969 vs 0.969, DNA 0.976 vs 0.976 | |
| &nbsp;&nbsp;`b_deploy` check | out-of-fold scores within 0.1 (max over proteins) | PASS | max difference 0.0000 | |
| &nbsp;&nbsp;`b_deploy` check | mean out-of-fold score difference within 0.01 | PASS | mean difference 0.0000 | |
| &nbsp;&nbsp;`b_deploy` check | applicability-domain cutoff within 0.01 | PASS | 0.8130 vs 0.8130 | |
<!-- VERIFY-GPU-TABLE-END -->

## 4. Not part of `make verify` or `make verify-gpu`

These need a download, an external install, or days of compute. Their inputs resolve from the archived data; the drivers are the `run_*.sh` files next to each stage.

| Result | Why not in `verify` | Entry point |
|---|---|---|
| New UniProt pull, substrate typing, UP clustering, structures, embeddings | network, GPU, hours | `make refresh` (`workflow/run_refresh.sh`), stages `pull` ... `embed_*` |
| MLP bake-off (16 splits), MTase-vs-UP training | GPU array / hours; the per-split results are archived and their statistics are verified | `2_training/run_bakeoff.sh`, `run_sam_vs_up.sh` |
| CLEAN and Horizyn runs | external installs (`workflow/setup_tools.sh`, `EXTERNAL_TOOLS.md`); the Horizyn scores are archived | `6_dark_cluster_annotation/run_clean_release.sh`, `run_horizyn_release.sh` |
| InterPro domain enrichment (Table S3) and the motif corpus | need the per-proteome UniProt tables with sequences (3 GB, not archived) | `7_motif_analysis/run_domain_table.sh`, `run_streme.sh`, `run_fimo_*.sh` |
| Figs 9 and 10 | ChimeraX | `figures/Fig9_10_structures/run_fig9_10.sh` |
| A new UMAP map | not reproducible run to run (same inputs, seed and environment give different points and 933-947 clusters), so `make verify` reuses the archived map (frozen); a new map: `UMAP_RECOMPUTE=1 bash 5_umap/run_umap.sh` | `5_umap/run_umap.sh` |

Run order for the stages that build on one another (the registry lists them in this order): clustering, benchmark reports, Fig 2, ensemble predictions, tree and Figs 4-5, cluster selection, dark-cluster tables and Fig 8, motifs, key numbers.
