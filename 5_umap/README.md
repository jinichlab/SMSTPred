# Stage 5 — Embedding landscape (UMAP) & cluster enrichment

Visualises the ProstT5 embedding space of all tree-of-life SAM-MTases and the
substrate-type / taxonomic structure within it, and zooms into family-specific clusters.

## Scripts

| Script | Purpose |
|---|---|
| `global_reumap_prostt5.py` | Global **UMAP** of every tree-of-life MTase from its ProstT5 embedding, coloured by predicted substrate type (the ProstT5 counterpart to the earlier ProtT5 map). Writes the 2-D coordinates, a static figure, and an interactive map. |
| `global_umap.py` | **Refresh release, ProtT5 (the map is frozen by default: UMAP and HDBSCAN are not reproducible from run to run, so an existing map of the same proteins is reused and only the scores and calls are refreshed; a release with other proteins gets a new map; `--recompute` / `UMAP_RECOMPUTE=1` forces one):** the same UMAP + HDBSCAN settings on the ProtT5 embeddings of a release's proteome proteins (all ensemble calls); writes `<release>/umap/` (coords CSV, static + interactive map, `umap_report.txt` with the neighbour purity). `RELEASE_DIR=<release> sbatch 5_umap/run_umap.sh`. |
| `cluster_selection.py` | **Refresh release:** the dark clusters (the definition of `dark_cluster_definition.py`), the metazoan protein-substrate and vertebrate DNA-substrate clusters (rules calibrated to reproduce the archived Fig 7A/B; `test_cluster_selection.py` checks that they recover the archived sets on the ProstT5 map), the per-cluster table and a crosswalk from the archived clusters to the new ones. `run_clusters.sh` also draws the Fig 7 grids (`figures/Fig7_clusters/build_cluster_grids.py --release`) into `<release>/umap/fig7/`. |
| `substrate_umap_datamap.py` | Substrate-type UMAP rendered with **datamapplot** — a clickable/interactive map (chosen over plotly), with SwissProt-anchored cluster labels. |
| `actino_metabolite_dark_clusters_prostt5.py` | Zoom-grid of the 14 Actinomycetota "dark-matter" clusters: zero SwissProt/reviewed members, ≥80% of members both Actinomycetota and predicted Metabolite (the rule is written out in `dark_cluster_definition.py`). |
| `plant_metabolite_dark_clusters_prostt5.py` | Zoom-grid of the 19 plant "dark-matter" clusters: zero SwissProt/reviewed members, ≥80% of members both plants (Viridiplantae) and predicted Metabolite. |

Both dark-cluster scripts hard-filter to the authoritative 33-cluster set in
`../6_dark_cluster_annotation/dark33_members.csv` (the same ground truth used by Tables S1 and S2
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
- `plant_metabolite_dark_clusters_prostt5.png` — the 19 plant small-molecule "dark-matter"
  clusters, **zero** SwissProt anchors by construction (see Stage 6 for CLEAN/Horizyn
  annotation of these).
- `actino_metabolite_dark_clusters_prostt5.png` — the 14 actinobacterial small-molecule
  "dark-matter" clusters, **zero** SwissProt anchors by construction (see Stage 6).
- `*_violin_prostt5.png` — per-branch substrate-load violins (vertebrate vs invertebrate
  DNA/Protein; plant-Metabolite; actinobacteria-Metabolite vs other bacteria). The two `*_other_violin_prostt5.png` files predate the class
  name Metabolite, still show the old label 'Other' in the picture, and have no generating script in this repository; the current
  versions of those panels are built by `figures/Fig4_trends/build_fig4_panels.py` (Fig 4C, 4D).

Substrate colours are consistent project-wide: RNA `#cb181d`, DNA `#6a51a3`,
Protein `#238b45`, Metabolite `#e6550d`.

## Why it matters

The map shows the embedding organises MTases into coherent, substrate- and family-specific
regions, recovers known families as labelled clusters, and surfaces **~thousands of
uncharacterised "dark-matter" enzymes** (plant + actinobacterial small-molecule MTases) as
distinct unanchored clusters — candidates for new biology, cross-validated against CLEAN
and Horizyn in Stage 6, and a motivation for the interpretability work in Stage 7.
