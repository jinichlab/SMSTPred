# Stage 5 — Embedding landscape (UMAP) & cluster enrichment

Visualises the ProstT5 embedding space of all tree-of-life SAM-MTases and the
substrate-type / taxonomic structure within it, and zooms into family-specific clusters.

## Scripts

| Script | Purpose |
|---|---|
| `global_reumap_prostt5.py` | Global **UMAP** of every tree-of-life MTase from its ProstT5 embedding, coloured by predicted substrate type (the ProstT5 counterpart to the earlier ProtT5 map). Writes the 2-D coordinates, a static figure, and an interactive map. |
| `substrate_umap_datamap.py` | Substrate-type UMAP rendered with **datamapplot** — a clickable/interactive map (chosen over plotly), with SwissProt-anchored cluster labels. |
| `actino_other_dark_clusters_prostt5.py` | Zoom-grid of the 14 Actinomycetota "dark-matter" clusters (Fig. 6D): zero SwissProt/reviewed members, ≥80% Actinomycetota & Other. |
| `plant_other_dark_clusters_prostt5.py` | Zoom-grid of the 19 plant "dark-matter" clusters (Fig. 6C): zero SwissProt/reviewed members, ≥80% plant & Other. |

Both dark-cluster scripts hard-filter to the authoritative 33-cluster set in
`../6_dark_cluster_annotation/dark33_members.csv` (the same ground truth used by Table 2
and Stage 6) and read `reviewed_status_full.csv` for reviewed/TrEMBL status and protein
names — **not** `uniprot_metadata.csv`, an earlier, incomplete lookup that undercounted
reviewed members and silently produced the wrong cluster membership (a stale 15-cluster
Actino figure and a stale 20-cluster plant figure). Both scripts assert on the authoritative
cluster set so this can't silently regress.

## Outputs

- `sam_mtase_umap_prostt5_coords.csv` — per-protein UMAP coordinates (+ labels).
- `sam_mtase_umap_prostt5_interactive.html` — interactive datamapplot.
- `sam_mtase_umap_prostt5_static_datamap.png` — static global map.
- `reviewed_status_full.csv` — per-Entry reviewed/TrEMBL status + protein name lookup
  (corrected; see caveat above). Input to the two dark-cluster scripts.

**Cluster-enrichment & lifestyle figures** (PNGs in this folder, built from the same UMAP /
HDBSCAN clustering and the deployment predictions):
- `metazoa_protein_clusters_prostt5.png` — metazoan protein-substrate clusters (PRMT/CARM),
  **anchored** on SwissProt members.
- `vertebrate_dna_clusters_prostt5.png` — vertebrate DNA-substrate clusters (DNMT3A/3B),
  **anchored** on SwissProt members.
- `plant_other_dark_clusters_prostt5.png` — the 19 plant small-molecule "dark-matter"
  clusters, **zero** SwissProt anchors by construction (see Stage 6 for CLEAN/Horizyn
  annotation of these).
- `actino_other_dark_clusters_prostt5.png` — the 14 actinobacterial small-molecule
  "dark-matter" clusters, **zero** SwissProt anchors by construction (see Stage 6).
- `*_violin_prostt5.png` — per-branch substrate-load violins (vertebrate vs invertebrate
  DNA/Protein; plant-Other; actinobacteria-Other vs other bacteria).

Substrate colours are consistent project-wide: RNA `#cb181d`, DNA `#6a51a3`,
Protein `#238b45`, Other `#e6550d`.

## Why it matters

The map shows the embedding organises MTases into coherent, substrate- and family-specific
regions, recovers known families as labelled clusters, and surfaces **~thousands of
uncharacterised "dark-matter" enzymes** (plant + actinobacterial small-molecule MTases) as
distinct unanchored clusters — candidates for new biology, cross-validated against CLEAN
and Horizyn in Stage 6, and a motivation for the interpretability work in Stage 7.
