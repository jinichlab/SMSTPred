# Stage 5 — Embedding landscape (UMAP) & cluster enrichment

Visualises the ProstT5 embedding space of all tree-of-life SAM-MTases and the
substrate-type / taxonomic structure within it, and zooms into family-specific clusters.

## Scripts

| Script | Purpose |
|---|---|
| `global_reumap_prostt5.py` | Global **UMAP** of every tree-of-life MTase from its ProstT5 embedding, coloured by predicted substrate type (the ProstT5 counterpart to the earlier ProtT5 map). Writes the 2-D coordinates, a static figure, and an interactive map. |
| `substrate_umap_datamap.py` | Substrate-type UMAP rendered with **datamapplot** — a clickable/interactive map (chosen over plotly), with SwissProt-anchored cluster labels. |

## Outputs

- `sam_mtase_umap_prostt5_coords.csv` — per-protein UMAP coordinates (+ labels).
- `sam_mtase_umap_prostt5_interactive.html` — interactive datamapplot.
- `sam_mtase_umap_prostt5_static_datamap.png` — static global map.

**Cluster-enrichment & lifestyle figures** (PNGs in this folder, built from the same UMAP /
HDBSCAN clustering and the deployment predictions, anchored on SwissProt members):
- `metazoa_protein_clusters_prostt5.png` — metazoan protein-substrate clusters (PRMT/CARM).
- `vertebrate_dna_clusters_prostt5.png` — vertebrate DNA-substrate clusters (DNMT3A/3B).
- `plant_other_clusters` / `plant_other_dark_clusters_prostt5.png` — plant small-molecule
  clusters (SABATH / phenylpropanoid) incl. **uncharacterised** ("dark") clusters with no
  SwissProt anchor.
- `actino_other_dark_clusters_prostt5.png` — actinobacterial small-molecule clusters,
  largely uncharacterised natural-product MTases.
- `*_violin_prostt5.png` — per-branch substrate-load violins (vertebrate vs invertebrate
  DNA/Protein; plant-Other; actinobacteria-Other vs other bacteria).

Substrate colours are consistent project-wide: RNA `#cb181d`, DNA `#6a51a3`,
Protein `#238b45`, Other `#e6550d`.

## Why it matters

The map shows the embedding organises MTases into coherent, substrate- and family-specific
regions, recovers known families as labelled clusters, and surfaces **~thousands of
uncharacterised "dark-matter" enzymes** (plant + actinobacterial small-molecule MTases) as
distinct unanchored clusters — candidates for new biology and a motivation for the
interpretability work in Stage 6.
