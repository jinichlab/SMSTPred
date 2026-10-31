# Reproducing the analyses

One entry point: the `Makefile`. Every step that runs on a CPU from the archived data is registered in
[`workflow/steps.py`](workflow/steps.py); `make verify` re-runs all of them in a scratch copy of the repository and compares what they
regenerate with the files committed here. The table in [section 3](#3-pass--fail-table) is the output of that check, from a fresh clone.

## 1. Quick start

```bash
git clone https://github.com/jinichlab/SMSTPred.git && cd SMSTPred
make env                          # conda environment "smstpred" from environment.yml (Python 3.12, pinned versions)
conda activate smstpred
make data DATA_SRC=zenodo:10.5281/zenodo.23655998      # download, extract and checksum the archived data into data/
#   already have the files?  make data DATA_SRC=/path/to/files LINK=--link
make check                        # SHA-256 of every file in data/ against MANIFEST.sha256
make verify                       # about 15 minutes on 16 CPU cores; writes verify/report.json and the table below
make verify-gpu                   # on a GPU machine: re-trains the models and compares outcomes (tolerances below)
```

| Target | What it does |
|---|---|
| `make env` | creates the conda environment (also the MEME suite, MMseqs2, HMMER and MAFFT from bioconda) |
| `make data` | puts the Zenodo files in `data/` (or `$SMSTPRED_DATA`) and extracts the archives; `DATA_SRC=zenodo:<doi>` or a local folder |
| `make check` | verifies the checksums of `data/` |
| `make tools` | downloads the Foldseek binary into `tools/` (only the bake-off baseline and SaProt embeddings need it) |
| `make steps` | lists every step: id, stage, folder, command |
| `make stage1` … `make stage7` | runs one stage in place (CPU steps) |
| `make figures` / `make tables` | runs the figure and table builders in place |
| `make verify` | runs every CPU step in a scratch copy, compares with the committed results, writes `verify/report.json` and this page's table |
| `make verify-gpu` | on a GPU machine: re-trains the deployment model, re-predicts all 1,259 proteomes, re-runs the MLP bake-off and the loss-curve sweeps, and compares the outcomes with stated tolerances (`workflow/gpu_checks.py`); uses `sbatch` for the loss-curve array when it exists. `sbatch workflow/slurm/verify_gpu.sh` does the same as one job |
| `make all` | `check`, every stage, `figures`, `verify` |

Scripts take their inputs from `data/` (or `$SMSTPRED_DATA`) and from other folders of this repository, and write their outputs
next to themselves. `make verify` never touches the committed files: it copies the repository (not `data/`, which is linked) to
`.verify_work/`, runs each step there, and compares every file the step wrote with the committed file of the same name.
A single step: `python workflow/run_steps.py --mode verify --only s2_table1`; one stage: `--stage 2`.

**The `run_*.sh` drivers** (GPU training, sweeps, Stage 7 searches) work on a SLURM cluster (`sbatch 2_training/run_train_deployment.sh`, from any
folder) and on any other machine (`bash 2_training/run_train_deployment.sh`; the `#SBATCH` lines are comments). Each starts in its own folder and
uses the conda environment that is already active, otherwise the `smstpred` environment from `make env` (`SMSTPRED_ENV=name` changes the name; the
job inherits the environment of the shell that submits it). Stage 6 is the exception: CLEAN and Horizyn are external installs with their own
environments (`6_dark_cluster_annotation/clean/environment_clean.yml`, `horizyn/environment_horizyn.yml`).

## 2. What "reproduces" means

| File type | Comparison | Tolerance |
|---|---|---|
| tables and statistics (`.csv`, `.tsv`, `.txt`, `.md`, `.json`) | same text outside the numbers; every number within tolerance | relative 1e-6 (1e-5 for the Stage 1 clustering, 1e-4 for the deployment-model predictions, which are float32 ensembles) |
| figures (`.png`) | same pixel size (within 2%); mean absolute grey-level difference | at most 12 of 255; this catches a missing panel, a wrong curve or a changed layout, not anti-aliasing or font hinting across matplotlib / FreeType builds |
| vector graphics, HTML, binary arrays and models | not compared | the numbers behind them are compared through the tables |

Statuses: **PASS** at least one regenerated file matched its committed version and none differed; **built** the step ran and wrote
its output but there is no committed counterpart to compare (figures ship as scripts, not images; the table of the manuscript
carries the strict checks); **FAIL** a step crashed, timed out, or a file differs beyond tolerance.

The strict, number-by-number checks are the tables behind the manuscript: Table 1, Table S4 (Holm-adjusted Wilcoxon),
the Friedman statistics of Fig 2A / Fig S4, the Stage 1 clustering tables, the dark-cluster definition (868 / 510 / 33), the
lineage comparisons, the prediction summary (112,019 proteins, 40,786 / 34,435 / 17,863 / 6,731 single calls, 38.2% outside
the applicability domain), the Horizyn validation (44.0% forward only, 85.7% either direction), Tables 2, S1, S2 and S3, and the
motif statistics (12,866 of 106,608 enzymes without a motif; 26 of 32 motifs at least two-fold class-specific).

## 3. Pass / fail tables

### CPU steps (`make verify`)

<!-- VERIFY-TABLE-START -->
Last verification: 2026-10-30T16:29:28 (Python 3.12.15, numpy 2.1.3, pandas 2.2.3, scipy 1.15.2, torch 2.8.0, matplotlib 3.10.9). 26 of 36 steps reproduce a committed result, 10 build their output without a committed reference to compare, 0 fail.

| Step | Reproduces | Result | Files compared | Time |
|---|---|---|---|---|
| `s1_clustering` | Stage 1 clustering tables at thresholds 0.5-0.95; Fig S1 | PASS | 4 image-match, 1 new-file, 2 numeric-match | 52 s |
| `s2_bakeoff_stats` | Friedman / Bonferroni-Dunn / Holm-Wilcoxon statistics of the benchmark (stats_report.txt) | PASS | 3 identical | 2 s |
| `s2_fig2a_zoom` | Fig 2A (cluster MCC, Friedman, Bonferroni-Dunn) and Fig S4 | PASS | 2 image-match | 5 s |
| `s2_sam_vs_up_stats` | Fig 2B statistics (Nadeau-Bengio paired test) | PASS | 2 identical | 3 s |
| `s2_sam_vs_up_plot` | Fig 2B bars | PASS | 1 identical | 2 s |
| `s2_fig2_main` | Fig 2 (A + B) | PASS | 1 identical | 3 s |
| `s2_table1` | Table 1 | PASS | 1 identical | 1 s |
| `s2_figs_S2_S5` | Figs S2 and S5 (MCC versus epoch) | PASS | 2 identical | 3 s |
| `s2_figs_S3_S6` | Figs S3 and S6 (train / validation loss) | PASS | 12 identical | 31 s |
| `s2_supp_stats` | Table S4 (Holm-adjusted Wilcoxon), Friedman tests for Fig 2A / Fig S4, the eleven ordinary iterations, per-iteration summary | PASS | 6 identical | 1 s |
| `s3_predict` | Deployment model: scores, calls, ensemble disagreement and applicability-domain flag for one organism | PASS | 1 identical, 1 new-file | 4 s |
| `s3_summary` | Prediction summary: 112,019 proteins, single / multi / no call, applicability domain, ensemble disagreement | PASS | 2 identical | 2 s |
| `s4_lineage_stats` | Lineage comparisons in the text (Micrococcales, Chlorophyta, plants, Actinomycetota; Mann-Whitney) | PASS | 2 identical | 1 s |
| `s5_fig6_umap` | Fig 6 (static UMAP) | PASS | 1 identical | 2 s |
| `s5_dark_definition` | Dark-cluster definition: 868 clusters, 510 without a reviewed member, 33 dark (3,478 proteins, 3.6%) | PASS | 2 identical | 1 s |
| `s6_dark_table` | Dark-cluster table (data behind Tables S1 and S2) | PASS | 1 identical | 0 s |
| `s6_horizyn_direction` | Horizyn validation: top-1 is a SAM-methylation reaction, 44.0% forward only vs 85.7% either direction (chance 2.3% vs 4.7%) | PASS | 2 identical | 10 s |
| `s7_domain_enrichment` | InterPro domain enrichment (Table S3 data) | PASS | 1 identical | 2 s |
| `s7_motif_table` | Table 2 (26 class-specific motifs) | PASS | 2 identical | 2 s |
| `s7_motif_status` | Per-protein motif status (12,866 of 106,608 motif-less) | PASS | 2 identical | 1 s |
| `s7_motif_specificity` | Motif class specificity (26 of 32 at least 2-fold) | PASS | 2 identical, 1 new-file | 4 s |
| `s7_family_reps` | 50%-identity family representatives (22,976) | PASS | 5 identical | 1 s |
| `s7_fig11_nomotif` | Fig 11 (motif-less enzymes on the UMAP) | PASS | 2 identical | 1 s |
| `f_fig3_widen` | Fig 3, step 1 (widen the iTOL gradient rings) | built (no reference) | 1 new-file | 1 s |
| `f_fig3_tree` | Fig 3, step 2 (tree + legend) | built (no reference) | 2 new-file | 8 s |
| `f_fig1` | Fig 1 (study overview) | built (no reference) | 3 new-file | 2 s |
| `f_fig4_panels` | Fig 4A-D | built (no reference) | 4 new-file | 3 s |
| `f_fig4_combined` | Fig 4 composite | built (no reference) | 1 new-file | 1 s |
| `f_fig5` | Fig 5 (proteome scaling) | built (no reference) | 1 new-file | 1 s |
| `f_fig7_grids` | Fig 7A-D cluster grids | built (no reference) | 4 new-file | 2 s |
| `f_fig7_combined` | Fig 7 composite | built (no reference) | 1 new-file | 1 s |
| `f_fig8_reaction_profiles` | Fig 8 (Horizyn reaction profiles of four dark clusters) | built (no reference) | 1 new-file | 4 s |
| `f_fig9_structures` | Structure panels (Figs 9 and 10) | built (no reference) | 2 new-file | 1 s |
| `f_table_S1` | Table S1 | PASS | 1 identical, 1 new-file | 1 s |
| `f_table_S2` | Table S2 | PASS | 1 identical, 1 new-file | 1 s |
| `f_table_S3` | Table S3 | PASS | 1 identical, 1 new-file | 1 s |
<!-- VERIFY-TABLE-END -->

### GPU re-runs (`make verify-gpu`)

A retrained network does not reproduce to the last digit across PyTorch / CUDA / GPU generations, so these steps compare outcomes with
tolerances set from measured noise. Training is deterministic within one environment (two runs of the deployment training gave identical weights), but the
committed models were trained in an earlier environment: re-training them with `environment.yml` on an RTX 6000 Ada moved the per-label thresholds by up to 0.008,
the out-of-fold scores by up to 0.012, the bake-off MCC by up to 0.04 per iteration (at most 0.005 per method), and **227 of 112,019 class calls (0.2%)**, so
that the re-trained model calls 40,904 proteins 'Metabolite' only (committed model: 40,786). The manuscript's numbers are those of the committed model in `2_training/deployment_model/`.

<!-- VERIFY-GPU-TABLE-START -->
Last verification: 2026-10-30T16:29:28 (Python 3.12.15, numpy 2.1.3, pandas 2.2.3, scipy 1.15.2, torch 2.8.0, matplotlib 3.10.9). 4 of 4 steps reproduce a committed result, 0 build their output without a committed reference to compare, 0 fail.

| Step | Reproduces | Result | Checks | Time |
|---|---|---|---|---|
| `b_deployment` | Deployment model re-trained from scratch: thresholds, out-of-fold MCC (0.969), out-of-fold scores | PASS | 4 of 4 passed | 13 s |
| `b_predict_all` | Predictions for all 1,259 proteomes (112,019 proteins) from the re-trained model | PASS | 4 of 4 passed | 32 s |
| `b_bakeoff_mlp` | Bake-off: the six MLP variants over the 16 cluster-disjoint iterations (Fig 2A, Table 1 rows) | PASS | 3 of 3 passed | 168 s |
| `b_loss_curves` | Train / validation loss curves behind Figs S3 and S6 (ten sweeps, SLURM array when sbatch exists) | PASS | 4 of 4 passed | 4206 s |
| &nbsp;&nbsp;`b_deployment` check | same proteins, clusters and folds | PASS | 7405 proteins | |
| &nbsp;&nbsp;`b_deployment` check | per-label thresholds within 0.02 | PASS | Metabolite 0.626 vs 0.634, Protein 0.635 vs 0.632, RNA 0.827 vs 0.826, DNA 0.109 vs 0.109 | |
| &nbsp;&nbsp;`b_deployment` check | per-label out-of-fold MCC within 0.01 | PASS | Metabolite 0.970 vs 0.970, Protein 0.971 vs 0.973, RNA 0.990 vs 0.990, DNA 0.945 vs 0.945 | |
| &nbsp;&nbsp;`b_deployment` check | out-of-fold scores within 0.03 | PASS | max difference 0.0117 | |
| &nbsp;&nbsp;`b_predict_all` check | same organism files | PASS | 1259 files | |
| &nbsp;&nbsp;`b_predict_all` check | identical class calls for at least 99.5% of proteins | PASS | 99.80% of 112,019 proteins | |
| &nbsp;&nbsp;`b_predict_all` check | scores within 0.05 | PASS | max difference 0.0228 | |
| &nbsp;&nbsp;`b_predict_all` check | applicability-domain flags identical | PASS |  | |
| &nbsp;&nbsp;`b_bakeoff_mlp` check | every MLP row of the re-run has a committed counterpart | PASS | 192 of 192 rows | |
| &nbsp;&nbsp;`b_bakeoff_mlp` check | max absolute MCC difference per iteration at most 0.06 | PASS | max 0.040, mean 0.0062 | |
| &nbsp;&nbsp;`b_bakeoff_mlp` check | mean MCC per method within 0.01 of the committed run | PASS | MLP_ESM-IF 0.9137 vs 0.9118; MLP_ProstT5 0.9563 vs 0.9514; MLP_ProstT5+ESM-IF 0.9483 vs 0.9474; MLP_ProtT5 0.9504 vs 0.9489; MLP_ProtT5+ESM-IF 0.9498 vs 0.9520; MLP_SaProt 0.9512 vs 0.9503 | |
| &nbsp;&nbsp;`b_loss_curves` check | a curves file for every committed sweep | PASS | 10 of 10 sweeps | |
| &nbsp;&nbsp;`b_loss_curves` check | median absolute validation-loss difference over all epochs at most 0.01 | PASS | median 0.0015, 95th percentile 0.1863 | |
| &nbsp;&nbsp;`b_loss_curves` check | median absolute final-epoch validation-loss difference at most 0.01 | PASS | median 0.0072, max 22.6662 | |
| &nbsp;&nbsp;`b_loss_curves` check | median absolute checkpoint-MCC difference at most 0.02 | PASS | median 0.0000, 95th percentile 0.0756 | |
<!-- VERIFY-GPU-TABLE-END -->

## 4. Not part of `make verify` or `make verify-gpu`

These steps need hours of compute, an external install, or a download that is not part of the Zenodo record. Their inputs
resolve from the archived data, and re-running them agrees with the committed results up to seed and hardware noise. The drivers are the
`run_*.sh` files next to each stage (see section 1 for how they pick their environment).

| Result | Why it is not in `verify` | Entry point |
|---|---|---|
| Embedding × loss × epoch sweeps | GPU, hours | `2_training/mlp_embedding_comparison.py`, `mlp_embedding_comparison_saprot.py`, `mlp_loss_comparison.py`, `mlp_loss_comparison_prostt5.py` |
| Classifier bake-off, baselines (profile HMM, MMseqs2, Foldseek, classical ML) | needs Foldseek (`make tools`), MMseqs2, HMMER; hours | `2_training/sam_classifier_comparison.py` (via `run_comparison.sh`); its six MLP rows are re-run by `make verify-gpu` (`run_bakeoff_mlp_only.sh`: same splits, the committed results were assembled that way, see `2_training/sam_classifier_comparison_output/PROVENANCE.md`). SAM vs UP: `sam_vs_up_prostt5_comparison.py` |
| Stage 1 UP clustering | the full 201,339-protein similarity matrix needs about 160 GB of disk | `1_clustering/up_clustering.py`; `clusters_UP_ProtT5_t0.5.csv` is provided instead |
| Stage 5 UMAP | long | `5_umap/global_reumap_prostt5.py` |
| Stage 4 tree, Stage 7 corpus | need the per-organism UniProt proteome tables (about 30 GB) | `4_phylo_trends/reference_tree.py`, `7_motif_analysis/1_build_corpus.py` |
| Stage 6 CLEAN / Horizyn runs | external installs | `CLEAN_APP_DIR`, `HORIZYN_DIR` (see `6_dark_cluster_annotation/README.md`); the Horizyn scores themselves are archived on Zenodo |
| Stage 7 motif search | MEME suite, MMseqs2 | `1_build_corpus.py` → `2_make_reps.py` → `run_streme.sh` → `run_fimo_all.sh` → `fimo_all_manifest.py` → `build_motif_table.py`; structure mapping: `run_fimo_cross.sh` → `fimo_cross_manifest.py` → `motif_structure_carriers.py` → `scan_pdb_ligands.py` (queries RCSB) → `structures/rebuild_from_carriers.py` |
| PyMOL snapshots (Figs 9, 10; Fig 1 structure panels) | PyMOL | the snapshots are included; the `view.pml` scripts in `7_motif_analysis/structures/` and `figures/Fig1_overview_schematic/render_prmt1_sam.pml` re-render them |

Run order for the figures that build on one another (the registry already lists them in this order): Fig 3 (`widen_rings.py`, then
`make_tree_figure.py`) before Fig 1, which embeds the tree panel; Fig 4 panels before the Fig 4 composite; Fig 7 grids before the Fig 7 composite.
