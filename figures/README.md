# Figure and table builders

Scripts that assemble the manuscript's composite figures and rendered tables. Figure
numbers follow the current manuscript draft. Analysis outputs they read live in the
numbered stage folders; scripts that sit in a stage folder instead are listed at the bottom.

Palette: the figure scripts use shades measured from recent Digital Discovery figures (blue `#4578aa`, green `#618d43`, orange `#ed8c44`, red `#c14445`, gold `#f4d23b`).
Class colours keep their hue: **RNA red `#c14445` · DNA purple `#6a51a3` · Protein green `#618d43` · Metabolite orange `#ed8c44`**. The UMAPs (Fig 6, Fig 11), the iTOL tree (Fig 3) and the
PyMOL structure panels (Figs 9 and 10) keep the original, more saturated shades (`#cb181d`, `#6a51a3`, `#238b45`, `#e6550d`) of the same hues. Fig 2 uses blue / green in panel A and orange / red in panel B
(never red next to green, which colour-blind readers cannot tell apart).

| Folder | Script(s) | Produces |
|---|---|---|
| `Fig1_overview_schematic/` | `build_fig1_flow.py`, `export_components.py` | Fig 1 pipeline overview (needs the Fig 3 `tree_only.png`; the PRMT1 structure panel is included in `biorender_components/`, re-render it with `render_prmt1_sam.pml` in PyMOL); `export_components.py` writes each sub-panel as a standalone transparent PNG for BioRender. |
| `Fig3_tree_landscape/` | `draw_tree.py` (new: draws the whole figure from the Newick + iTOL files, no iTOL), `make_tree_figure.py`, `widen_rings.py` | Fig 3 circular tree from the iTOL SVG export. `widen_rings.py` rewrites the five gradient-ring wedge paths in the SVG to a wider band (default 3x); `make_tree_figure.py` rasterizes, crops, adds the legend and clade highlights. |
| `Fig4_trends/` | `build_fig4_panels.py`, `make_fig4_trends.py` | Fig 4A–D violin panels from the per-organism predictions and organism metadata (`build_fig4_panels.py`), then the 2×2 composite (`make_fig4_trends.py`). Computes the Kruskal-Wallis / Mann-Whitney statistics quoted in the text. For a refresh release: `RELEASE_DIR=<release> sbatch figures/Fig4_trends/run_fig4.sh` (`build_fig4_panels.py --release`, then `make_fig4_trends.py --dir`), writing `<release>/fig4/` with `fig4_statistics.csv` (every test) and `fig4_organism_table.csv`. |
| `Fig5_proteome_scaling/` | `scaling_by_clade.py` (Fig 5: eukaryotes fitted clade by clade, with a model comparison), `scaling_loglog.py` (lineage-vs-rest fits, kept as the previous version), `proteome_vs_sam_domain_highlight.py` | Fig 5 SAM-MTase count vs proteome size with per-lineage OLS fits and the fold-faster-scaling ratios. Reads `proteome_vs_sam_bylabel.csv` (one row per organism: domain, proteome size, SAM-MTase counts total and per substrate label; built from the 1,259 per-organism prediction CSVs in `3_predictions/predictions_output/` plus reference-proteome sizes). The archived table is in this folder; for a refresh release `build_proteome_vs_sam_table.py` builds it from `<release>/proteome_summary.csv` + `predictions_all_proteomes.csv.gz` and `run_fig5.sh` (`RELEASE_DIR=<release> sbatch figures/Fig5_proteome_scaling/run_fig5.sh`) draws the figure into `<release>/fig5/` with `fig5_fits.csv` (slopes, R2, ratio with a bootstrap 95% interval). |
| `Fig7_clusters/` | `build_cluster_grids.py`, `make_fig7_clusters.py` | Fig 7A–D cluster-UMAP grids for the four dark clusters c110/c551/c667/c714 (from `5_umap/sam_mtase_umap_prostt5_coords.csv`, each outlined in its own colour) and the composite. |
| `Fig8_reaction_profiles/` | `make_fig8_reaction_profiles.py` | Fig 8: for each of the four clusters, the cluster-mean Horizyn score of every Rhea reaction (direction-collapsed), the SAM-methylation reactions marked, and one reaction drawn as substrate → product (the one most members voted for; for c714 the furocoumarin reaction CLEAN names, which is the second most voted). Reads `horizyn_dark_scores.npz` from `data/`. With `--release <release>` (the refresh chain) it draws the clusters chosen by `select_panels.py` (CLEAN reliable + methyltransferase, and CLEAN's EC inside Horizyn's dominant reaction family; on release 2026_03 c317, c353, c512, c729), reads the release's Horizyn scores, takes the panel titles from `panel_titles.csv` and writes `<release>/fig8/` (also `fig8_clusters.csv`, which Fig 7 outlines). Rhea entries are mapped to their master reaction with `rhea-directions.tsv`: the master is not a fixed offset from its entries. |
| `Fig9_10_structures/` | `make_panels.py` | Lays the PyMOL motif-on-structure snapshots (`Known/`, `Distal/`; rendered from the `view.pml` scripts in `7_motif_analysis/structures/`) out as labelled A–D panels (Fig 9 known sites, Fig 10 distal sites). For a release (`run_fig9_10.sh`, part of the refresh chain): `select_structures.py` picks four motifs for each figure by rule from `motif_placements.csv` (Fig 9 = the closest SAM + substrate motif of each class, Fig 10 = the furthest distal motif of each class; a class without one gets the next closest / furthest of any class; weak picks are skipped), `render_panels.py` renders them headlessly with ChimeraX, `make_panels.py --release` lays them out. |
| `FigS1_similarity_thresholds/` | `make_figS1.py` | Combines the two cosine-similarity histograms written by `1_clustering/sam_clustering.py`. |
| `TableS1_actinomycetota_dark_clusters/` | `build_table.py` | Table S1 rendered from `6_dark_cluster_annotation/table2_dark_clusters.csv`. |
| `TableS2_plant_dark_clusters/` | `build_table.py` | Table S2, same source, plant clusters. |
| `TableS3_domain_enrichment/` | `build_domain_enrichment_table.py` | Table S3 from `7_motif_analysis/domain_enrichment.csv`. |

## Figure scripts that live in a stage folder

| Figure / table | Script |
|---|---|
| Fig 2, Fig S4, Table S4 | `2_training/plot_fig2_main.py`, `plot_cluster_mcc_zoom.py`, `classifier_mcc_stats.py`, `sam_vs_up_prostt5_stats.py`, `plot_sam_vs_up_prostt5.py`, `benchmark_supplement_stats.py` |
| Table 1 | `2_training/build_metrics_table.py` |
| Figs S2, S5 | `2_training/plot_mcc_curves_grid.py` (data: `mlp_embedding_comparison.py`, `mlp_embedding_comparison_saprot.py`, `mlp_loss_comparison_prostt5.py`) |
| Figs S3, S6 | `2_training/plot_train_val_loss_curves.py` (data: `train_val_loss_curves.py`; raw curves on Zenodo) |
| Fig 4A (original) | `4_phylo_trends/protein_cellularity_violin.py` |
| Fig 6 | `5_umap/plot_umap_static.py` |
| Table 2, Fig 11 | `7_motif_analysis/build_motif_table.py`, `no_motif_umap.py` |

## Notes

- Small input tables are in [`../data/`](../data); the large embeddings are archived on Zenodo
  (see [`../data/README.md`](../data/README.md)); scripts find them relative to the repository or via
  `$SMSTPRED_DATA`.
- `make verify` (see [`../REPRODUCING.md`](../REPRODUCING.md)) reruns every builder here from a fresh clone and compares the
  tables with the committed versions (numbers within a stated tolerance). Figures are checked for the same size and content
  (mean grey-level difference of at most 12/255); pixels can differ slightly from the manuscript PNGs (anti-aliasing, label
  offsets of a pixel) across matplotlib/FreeType builds. Table S1's line wraps need Python 3.12 for an exact match (under 3.9,
  `textwrap` breaks one long chemical name at different points).
- The Fig 3 raw iTOL SVG export (`V4prbvTzDgPu09RO61xuyg.svg`) is an input to
  `Fig3_tree_landscape/` and is archived on Zenodo (as `itol_tree_export_Fig3.svg`), not in git.
