# Reproducing the analyses

Status of every step, from a fresh clone. "Verified" means it was run from a clean `git clone` in an
environment built from `environment.yml`, using only this repository and the Zenodo data (see
[`data/README.md`](data/README.md)).

## Setup

```bash
git clone https://github.com/jinichlab/SMSTPred.git && cd SMSTPred
conda env create -f environment.yml && conda activate classifier_env   # Python 3.12, pinned versions

# Zenodo data -> data/ (or any folder, then: export SMSTPRED_DATA=/path/to/folder)
cd data && for t in embeddings_tree_of_life_*.tar; do tar -xf "$t"; done && cd ..
```

Scripts take their inputs from `data/` (or `$SMSTPRED_DATA`) and from other folders of this repository, and write
their outputs next to themselves. The environment was validated with the CPU build of PyTorch 2.8; GPU work
needs a CUDA-enabled install. Not in the environment: the MEME suite (STREME/FIMO), MMseqs2, Foldseek, HMMER/MAFFT
(bake-off baselines), and the CLEAN and Horizyn checkouts (Stage 6).

## A. Verified: figures, tables and predictions from the archived data

Run each from the script's own folder, e.g. `cd figures/Fig4_trends && python build_fig4_panels.py`.

| Result | Script(s) |
|---|---|
| Stage 3 predictions (any ProstT5 embedding file) | `3_predictions/predict.py --input <embeddings.h5> --out-dir <dir>` |
| Fig 2, Fig S4 and their statistics | `2_training/classifier_mcc_stats.py`, `plot_cluster_mcc_zoom.py`, `sam_vs_up_prostt5_stats.py`, `plot_sam_vs_up_prostt5.py`, `plot_fig2_main.py` |
| Table 1 | `2_training/build_metrics_table.py` |
| Figs S2, S3, S5 | `2_training/plot_mcc_curves_grid.py` |
| Fig 3 (two steps, in this order) | `figures/Fig3_tree_landscape/widen_rings.py`, then `make_tree_figure.py` |
| Fig 4A–D and composite | `figures/Fig4_trends/build_fig4_panels.py`, then `make_fig4_trends.py` (4A original: `4_phylo_trends/protein_cellularity_violin.py`) |
| Fig 5 | `figures/Fig5_proteome_scaling/proteome_vs_sam_domain_highlight.py` |
| Fig 6 | `5_umap/plot_umap_static.py` |
| Fig 7A–E and composite | `figures/Fig7_clusters/build_cluster_grids.py`, `render_reaction_structures.py`, then `make_fig7_clusters.py` (dark-cluster zoom grids: `5_umap/plant_other_…`, `actino_other_…`) |
| Fig 10 (static and interactive) | `7_motif_analysis/no_motif_umap.py`, `no_motif_umap_interactive.py` |
| Tables S1, S2 (from the Stage 6 table, rebuilt by `6_dark_cluster_annotation/build_table2.py`) | `figures/TableS1_…/build_table.py`, `figures/TableS2_…/build_table.py` |
| Table S3 | `figures/TableS3_domain_enrichment/build_domain_enrichment_table.py` |

**What "reproduces" means here.** Every table CSV and statistics file regenerated identical to the committed
versions. Figures regenerated with the same content (same curves, values and labels); pixels differ slightly from
the manuscript PNGs (anti-aliasing, label offsets of a pixel) across matplotlib/FreeType builds.

## B. Not reproducible from the repository and Zenodo alone

| Result | Why not | What would be needed |
|---|---|---|
| Table 2 (`7_motif_analysis/build_motif_table.py`) | reads `fimo_all_*/fimo.tsv` and `class_*.fasta` | the Stage 7 intermediates (FIMO scan of all class members, class FASTAs) |
| Fig 1 (`figures/Fig1_overview_schematic/`) | composes PyMOL renders and BioRender components that are not shipped (the tree panel comes from the Fig 3 step) | the two structure images |
| Fig S1B (deployment-domain similarity, `1_clustering/sam_clustering.py` validation stage) | needs `prediction_files_mlp/` | not archived (listed in `data/README.md`) |
| UP clustering (`1_clustering/up_clustering.py`) | needs `prediction_files_mlp/` and the UP ProtT5 embeddings (`per-protein.h5`) | not archived (listed in `data/README.md`) |
| Stage 4 tree (`4_phylo_trends/reference_tree.py`), Stage 7 corpus (`7_motif_analysis/1_build_corpus.py`) | need the per-organism UniProt proteome tables (about 30 GB) | download of the 1,259 UniProt reference proteomes |
| Foldseek baseline (part of `2_training/sam_full_vs_balanced_classifier_comparison.py`) | needs AlphaFold models of the 7,405 SAM proteins (`sam_structures/`) and the `foldseek` binary | AlphaFold DB download |
| Stage 6 CLEAN / Horizyn runs | external installs | `CLEAN_APP_DIR`, `HORIZYN_DIR` (see `6_dark_cluster_annotation/README.md`) |

## C. Heavy computations (path-checked, not re-run here)

These need a GPU and hours, and were not re-run for this check; their inputs resolve from the archived data and
re-training will agree with the committed results up to seed and hardware noise. SLURM drivers are the
`run_*.sh` files next to each stage.

- Embedding x loss x epoch sweeps: `2_training/mlp_embedding_comparison.py`, `mlp_embedding_comparison_saprot.py`,
  `mlp_loss_comparison.py`, `mlp_loss_comparison_prostt5.py`
- Classifier bake-off: `2_training/sam_full_vs_balanced_classifier_comparison.py`; SAM vs UP:
  `sam_vs_up_prostt5_comparison.py`
- Deployment model: `2_training/train_deployment_ensemble.py` (rewrites `2_training/deployment_model/`)
- Stage 1 clustering: `1_clustering/sam_clustering.py` (the deployment-domain validation stage is skipped when
  `prediction_files_mlp/` is absent)
- Stage 5 UMAP: `5_umap/global_reumap_prostt5.py`
- Stage 7 motif pipeline, in order: `1_build_corpus.py` -> `2_make_reps.py` (MMseqs2) -> `run_streme.sh` (STREME)
  -> `run_fimo_all.sh` -> `fimo_all_manifest.py` -> `build_motif_table.py`; structure mapping:
  `run_fimo_cross.sh` -> `fimo_cross_manifest.py` -> `motif_structure_carriers.py` -> `scan_pdb_ligands.py`
  (queries RCSB) -> `structures/rebuild_from_carriers.py`
