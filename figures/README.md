# Figure and table builders

Scripts that assemble the manuscript's composite figures and rendered tables. Figure
numbers follow the current manuscript draft. Analysis outputs they read live in the
numbered stage folders; scripts that sit in a stage folder instead are listed at the bottom.

Palette throughout: **RNA `#cb181d` · DNA `#6a51a3` · Protein `#238b45` · Other `#e6550d`**.

| Folder | Script(s) | Produces |
|---|---|---|
| `Fig1_overview_schematic/` | `build_fig1_flow.py`, `export_components.py` | Fig 1 pipeline overview; `export_components.py` writes each sub-panel as a standalone transparent PNG for BioRender. |
| `Fig3_tree_landscape/` | `make_tree_figure.py`, `widen_rings.py` | Fig 3 circular tree from the iTOL SVG export. `widen_rings.py` rewrites the five gradient-ring wedge paths in the SVG to a wider band (default 3x); `make_tree_figure.py` rasterizes, crops, adds the legend and clade highlights. |
| `Fig4_trends/` | `build_fig4_panels.py`, `make_fig4_trends.py` | Fig 4A–D violin panels from the per-organism predictions and organism metadata (`build_fig4_panels.py`), then the 2×2 composite (`make_fig4_trends.py`). Computes the Kruskal-Wallis / Mann-Whitney statistics quoted in the text. |
| `Fig5_proteome_scaling/` | `proteome_vs_sam_domain_highlight.py` | Fig 5 SAM-MTase count vs proteome size with per-lineage OLS fits and the fold-faster-scaling ratios. Reads `proteome_vs_sam_bylabel.csv` (one row per organism: domain, proteome size, SAM-MTase counts total and per substrate label; built from the 1,259 per-organism prediction CSVs in `3_predictions/predictions_output/` plus reference-proteome sizes). The table itself is included in this folder; the script that builds it is not yet in the repository. |
| `Fig7_clusters/` | `build_cluster_grids.py`, `render_reaction_structures.py`, `make_fig7_clusters.py` | Fig 7A–D cluster-UMAP grids (from `5_umap/sam_mtase_umap_prostt5_coords.csv`), Fig 7E reaction schemes for clusters c551/c667/c110/c714 (RDKit), and the composite. |
| `Fig9_structures/` | `make_panels.py` | Lays the PyMOL motif-on-structure snapshots (`7_motif_analysis/structures/`) out as labelled A–D panels (Figs 8–9). |
| `FigS1_similarity_thresholds/` | `make_figS1.py` | Combines the two cosine-similarity histograms written by `1_clustering/sam_clustering.py`. |
| `TableS1_actinomycetota_dark_clusters/` | `build_table.py` | Table S1 rendered from `6_dark_cluster_annotation/table2_dark_clusters.csv`. |
| `TableS2_plant_dark_clusters/` | `build_table.py` | Table S2, same source, plant clusters. |
| `TableS3_domain_enrichment/` | `build_domain_enrichment_table.py` | Table S3 from `7_motif_analysis/domain_enrichment.csv`. |

## Figure scripts that live in a stage folder

| Figure / table | Script |
|---|---|
| Fig 2, Fig S4 | `2_training/plot_fig2_main.py`, `plot_cluster_mcc_zoom.py`, `classifier_mcc_stats.py`, `sam_vs_up_prostt5_stats.py`, `plot_sam_vs_up_prostt5.py` |
| Table 1 | `2_training/build_metrics_table.py` |
| Figs S2, S3, S5 | `2_training/plot_mcc_curves_grid.py` (data: `mlp_embedding_comparison.py`, `mlp_embedding_comparison_saprot.py`, `mlp_loss_comparison_prostt5.py`) |
| Fig 4A (original) | `4_phylo_trends/protein_cellularity_violin.py` |
| Fig 6 | `5_umap/plot_umap_static.py` |
| Table 2, Fig 10 | `7_motif_analysis/build_motif_table.py`, `no_motif_umap.py` |

## Notes

- Small input tables are in [`../data/`](../data); the large embeddings are archived on Zenodo
  (see [`../data/README.md`](../data/README.md)); scripts find them relative to the repository or via
  `$SMSTPRED_DATA`.
- Checked against the manuscript versions, these reproduce pixel for pixel from the repository: Fig 4A–D,
  Fig 7A–D and Tables 1, S1, S2, S3. Table S1 needs Python 3.12 (matplotlib 3.10) for an exact match: under Python
  3.9, `textwrap` breaks one long chemical name at different points, so only that table's line wraps differ.
- The Fig 3 raw iTOL SVG export (`V4prbvTzDgPu09RO61xuyg.svg`) is an input to
  `Fig3_tree_landscape/` and is archived on Zenodo (as `itol_tree_export_Fig3.svg`), not in git.
