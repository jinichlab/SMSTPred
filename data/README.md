# Data

Small tables the analysis needs, kept in git. The large inputs (precomputed embeddings, the raw tree
export) are archived on **Zenodo** (doi:10.5281/zenodo.23655998; reserved, resolves once the record is published); see the record's README for layouts,
checksums and provenance.

## In this folder

| File | Rows | What it is |
|---|---|---|
| `labels_SAM_full.csv` | 7,405 | Curated SAM-methyltransferase training set (UniProtKB/Swiss-Prot, reviewed): `Entry, Reviewed, Entry Name, Protein names, Length, Sequence, Annotation, InterPro, Catalytic activity, RheaID, type, ChEBI, Metabolite, Protein, RNA, DNA`. Label counts: Metabolite 1,290 · Protein 1,015 · RNA 5,006 · DNA 108 (multi-label, so counts are positive labels, not a partition). |
| `labels_SAM_balanced.csv` | 2,924 | Class-balanced subset of the above: `Entry, Sequence, Annotation, RNA, DNA, Metabolite, Protein`. Metabolite 689 · Protein 1,015 · RNA 1,126 · DNA 108. Not used in the manuscript; the optional balanced-training arm of `2_training/sam_classifier_comparison.py` reads it, and the classifier comparison reads it in every run to keep its random splits unchanged. |
| `labels_UP.csv` | 201,339 | Broad enzyme comparison set (UP): `Entry` plus the four labels only; sequences are retrievable from UniProt by accession. Metabolite 133,297 · Protein 30,967 · RNA 36,846 · DNA 5,691. |
| `selected_1259_organisms.csv` | 1,259 | The tree-of-life reference proteomes: proteome id, organism, taxonomic lineage, BUSCO score, lifestyle and its confidence, genus, and the metazoa / vertebrate / parasite flags used by the Stage 4 and 5 trend analyses. |

Column meaning for the labels: `Metabolite` = a substrate that is not a protein, RNA or DNA: small molecules and metabolites, including lipids and
polyprenyl quinones (the UP set also holds some oligo- and polysaccharides), `Protein`, `RNA`, `DNA`; 1 means
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
AlphaFold DB models plus 26 locally predicted ones), the raw iTOL SVG used for Fig 3, the raw train/validation loss
curves (`train_val_loss_curves.tar`), and the Horizyn scores behind Fig 8 (`horizyn_dark_scores.npz`,
`horizyn_dark_top50.csv`, `horizyn_rhea_reaction_embeddings.npz`).

### Provenance

| Item | Source |
|---|---|
| SAM training set (`labels_SAM_full.csv`) | UniProtKB/Swiss-Prot release **2024_05** |
| UP comparison set (`labels_UP.csv`) | UniProtKB release **2024_04** |
| Tree-of-life reference proteomes (1,259 organisms) | UniProt release **2026_01**; the per-organism tables of 146 organisms were re-pulled on 2026-09-08 (release 2026_03) to refresh the Reviewed flag. How the organisms were chosen is in the `selection_reason` column of `selected_1259_organisms.csv` |
| ProtT5 embeddings of the SAM and UP sets | `Rostlab/prot_t5_xl_half_uniref50-enc`, mean over residues, stored float16 |
| ProtT5 embeddings of the tree-of-life set | `Rostlab/prot_t5_xl_uniref50`, float32, masked mean over tokens (the 608 proteins that occur in both pipelines agree to cosine ≥ 0.9995) |
| ProstT5 embeddings (all sets) | `Rostlab/ProstT5`, `<AA2fold>` prefix, fp16, mean over residue tokens, stored float16 |
| ESM-IF embeddings | ESM-IF1 (AlphaFold DB weights, v6), chain A, 512-D, fair-esm 2.0.0 |
| SaProt embeddings | `westlake-repl/SaProt_650M_AF2` with Foldseek 3Di strings |
| Structures | AlphaFold DB models; the 26 proteins without one were predicted with ColabFold 1.6.1 (AlphaFold2-ptm, five models, MSAs from the ColabFold MMseqs2 server, no templates or relaxation; the highest-pLDDT model was kept) |
| Protein length | sequences of at most 1,000 residues were embedded |
| Tree-of-life protein counts | ProstT5: 112,019 (the final prediction set). ProtT5 archive: 111,343. `prediction_files_mlp.tar` (earlier ProtT5-MLP run, used for Fig S1B): 111,337, all among the 112,019; 608 of them are training proteins and have similarity 1.0 to themselves |

## Using the archived data

Every script locates its inputs relative to the repository, so a fresh clone runs without edits.
Large inputs are looked up in one data directory, **`data/` by default**, or wherever the environment
variable `SMSTPRED_DATA` points:

```bash
make data DATA_SRC=zenodo:10.5281/zenodo.23655998     # download, extract and checksum
make data DATA_SRC=/path/to/files LINK=--link         # or: use a local copy of the record (symlinked)
make check                                            # SHA-256 of every file against MANIFEST.sha256
export SMSTPRED_DATA=/elsewhere                       # only if the data should live outside <repo>/data
```

`make data` extracts the archives where the scripts expect them: the tree-of-life embeddings, `prediction_files_mlp/` and
`sam_structures/` into `data/`, `motif_analysis_intermediates.tar` into `7_motif_analysis/`, and `train_val_loss_curves.tar` into
`2_training/`.

Files are expected under their Zenodo names (`labels_SAM_full.csv`, `embeddings_SAM_ProstT5.h5`, ...).
Paths to other repository folders (predictions, UMAP coordinates, result tables) are derived from each
script's own location.

## Input that is not archived

Needed only to re-run the two steps listed; neither is required for the figures and tables built from the files
above.

| Input (expected under `SMSTPRED_DATA`) | Needed by | Where it comes from |
|---|---|---|
| `reference_proteomes/` (per-organism UniProt proteome tables, about 30 GB) | `4_phylo_trends/reference_tree.py`, `7_motif_analysis/1_build_corpus.py` | UniProt reference proteomes of the organisms in `selected_1259_organisms.csv` |

`reference_proteomes/` can be rebuilt (from the current UniProt, not the release the paper used) with `make pull`
([`0_uniprot_pull/`](../0_uniprot_pull)); the per-organism files it writes have the same names and columns.

External tools and installs: `foldseek` and `mmseqs` on `PATH` (override with `FOLDSEEK_BIN`); for Stage 6, a
CLEAN checkout (`CLEAN_APP_DIR`, used by `6_dark_cluster_annotation/clean/run_clean_dark33*.sh`) and a
Horizyn checkout (`HORIZYN_DIR`, used by the Horizyn scripts). The Rhea reaction tables Horizyn reads are
included in `6_dark_cluster_annotation/horizyn/inputs/` (override with `SMSTPRED_HORIZYN_INPUTS`).

## License

The tables in this folder are licensed CC BY 4.0, like the Zenodo record. They derive from UniProtKB (with Rhea and ChEBI
annotation), which is CC BY 4.0: please credit UniProt when reusing them.
