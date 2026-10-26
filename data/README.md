# Data

Small tables the analysis needs, kept in git. The large inputs (precomputed embeddings, the raw tree
export) are archived on **Zenodo** (DOI to be added at publication); see the record's README for layouts,
checksums and provenance.

## In this folder

| File | Rows | What it is |
|---|---|---|
| `labels_SAM_full.csv` | 7,405 | Curated SAM-methyltransferase training set (UniProtKB/Swiss-Prot, reviewed): `Entry, Reviewed, Entry Name, Protein names, Length, Sequence, Annotation, InterPro, Catalytic activity, RheaID, type, ChEBI, Other, Protein, RNA, DNA`. Label counts: Other 1,290 · Protein 1,015 · RNA 5,006 · DNA 108 (multi-label, so counts are positive labels, not a partition). |
| `labels_SAM_balanced.csv` | 2,924 | Class-balanced subset of the above (Fig S3): `Entry, Sequence, Annotation, RNA, DNA, Other, Protein`. Other 689 · Protein 1,015 · RNA 1,126 · DNA 108. |
| `labels_UP.csv` | 201,339 | Broad enzyme comparison set (UP): `Entry` plus the four labels only; sequences are retrievable from UniProt by accession. Other 133,297 · Protein 30,967 · RNA 36,846 · DNA 5,691. |
| `selected_1259_organisms.csv` | 1,259 | The tree-of-life reference proteomes: proteome id, organism, taxonomic lineage, BUSCO score, lifestyle and its confidence, genus, and the metazoa / vertebrate / parasite flags used by the Stage 4 and 5 trend analyses. |

Column meaning for the labels: `Other` = small-molecule / metabolite substrates, `Protein`, `RNA`, `DNA`; 1 means
the enzyme acts on that substrate class, assigned from curated catalytic-activity annotation (reaction, Rhea,
ChEBI). See the manuscript Methods.

## Data that lives elsewhere in the repository

| Data | Location |
|---|---|
| Per-organism substrate predictions for all 1,259 proteomes (about 112,000 methyltransferases) | `3_predictions/predictions_output/` |
| Deployment model weights and frozen thresholds | `2_training/deployment_model/` |
| UMAP coordinates, interactive map | `5_umap/` |
| SAM clusterings at cosine thresholds 0.5–0.95 | `1_clustering/sam_clustering_output/` |
| Fig 5 input (per-organism proteome size and SAM-methyltransferase counts) | `figures/Fig5_proteome_scaling/proteome_vs_sam_bylabel.csv` |
| CLEAN / Horizyn results, dark-cluster table | `6_dark_cluster_annotation/` |
| STREME / FIMO motifs, domain enrichment, structure placements | `7_motif_analysis/` |

## On Zenodo (not in git)

Embeddings (ProtT5, ProstT5, ESM-IF, SaProt for the SAM set; ProstT5 for UP; ProstT5 and ProtT5 for the
tree-of-life set), the UP clustering (`clusters_UP_ProtT5_t0.5.csv`), copies of the label tables above, the
per-organism ProtT5-MLP prediction files (`prediction_files_mlp.tar`), the Stage 7 intermediates
(`motif_analysis_intermediates.tar`: single-label corpus, class and representative FASTAs, FIMO hit tables), the UP
ProtT5 embeddings (`embeddings_UP_ProtT5.h5`), the structures of the SAM proteins (`sam_structures.tar.gz`: 7,379
AlphaFold DB models plus 26 locally predicted ones), and the raw iTOL SVG used for Fig 3.

## Using the archived data

Every script locates its inputs relative to the repository, so a fresh clone runs without edits.
Large inputs are looked up in one data directory, **`data/` by default**, or wherever the environment
variable `SMSTPRED_DATA` points:

```bash
# download the Zenodo record into data/ (or any folder) and extract the archives
cd data
tar -xf embeddings_tree_of_life_ProstT5.tar     # -> tree_of_life_ProstT5/
tar -xf embeddings_tree_of_life_ProtT5.tar      # -> tree_of_life_ProtT5/
tar -xf prediction_files_mlp.tar                # -> prediction_files_mlp/
tar -xzf sam_structures.tar.gz                  # -> sam_structures/ (Foldseek baseline)
export SMSTPRED_DATA=$PWD                       # only needed if the data live outside <repo>/data
cd .. && tar -xf data/motif_analysis_intermediates.tar   # Stage 7 corpus, FASTAs, FIMO tables -> 7_motif_analysis/
```

Files are expected under their Zenodo names (`labels_SAM_full.csv`, `embeddings_SAM_ProstT5.h5`, ...).
Paths to other repository folders (predictions, UMAP coordinates, result tables) are derived from each
script's own location.

## Input that is not archived

Needed only to re-run the two steps listed; neither is required for the figures and tables built from the files
above.

| Input (expected under `SMSTPRED_DATA`) | Needed by | Where it comes from |
|---|---|---|
| `reference_proteomes/` (per-organism UniProt proteome tables, about 30 GB) | `4_phylo_trends/reference_tree.py`, `7_motif_analysis/1_build_corpus.py` | UniProt reference proteomes of the organisms in `selected_1259_organisms.csv` |

External tools and installs: `foldseek` and `mmseqs` on `PATH` (override with `FOLDSEEK_BIN`); for Stage 6, a
CLEAN checkout (`CLEAN_APP_DIR`, used by `6_dark_cluster_annotation/clean/run_clean_dark33*.sh`) and a
Horizyn checkout (`HORIZYN_DIR`, used by the Horizyn scripts). The Rhea reaction tables Horizyn reads are
included in `6_dark_cluster_annotation/horizyn/inputs/` (override with `SMSTPRED_HORIZYN_INPUTS`).
