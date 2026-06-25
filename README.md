# Predicting methyltransferase substrate specificity across the tree of life with structure-aware protein language models

Analysis code and results for the manuscript. The project builds a classifier that
predicts the **substrate type** of a SAM-dependent methyltransferase (MTase) —
**DNA / RNA / Protein / Other (small-molecule)** — directly from a protein's sequence via
a **structure-aware protein language model (ProstT5)**, validates it against
homology and classical-ML baselines, deploys it across the tree of life, and interprets
*what* the model reads to make each call.

## The pipeline (stages)

| Stage | Folder | What it does |
|---|---|---|
| 1 | [`1_clustering/`](1_clustering) | Homology-aware **cluster-disjoint** partitioning of the training sets (greedy cosine clustering + threshold selection) so cross-validation never leaks homologues. |
| 2 | [`2_training/`](2_training) | Loss/epoch tuning, the full **classifier bake-off** (embeddings vs homology baselines vs classical ML), the **SAM-vs-UP** training-set comparison, and the final **deployment model**. |
| 3 | [`3_predictions/`](3_predictions) | Scores new proteins with the deployment model (per-label score, call, ensemble disagreement, applicability-domain flag). |
| 4 | [`4_phylo_trends/`](4_phylo_trends) | Phylogenetic / taxonomic trends in the predicted substrate usage across the tree of life. |
| 5 | [`5_umap/`](5_umap) | Global ProstT5 UMAP of all tree-of-life MTases + substrate-coloured maps and cluster-enrichment figures. |
| 6 | [`6_motif_analysis/`](6_motif_analysis) | **Interpretability**: discriminative sequence motifs (STREME), per-residue attribution (Integrated Gradients), and structural mapping of both — *why* the classifier works. |

Each folder has its own `README.md` with scripts, inputs, outputs, and parameters.

## Headline results

- **ProstT5 is the best embedding.** Under cluster-disjoint CV the ProstT5 MLP reaches
  **MCC 0.955**, beating sequence-homology (MMseqs2 1-NN 0.917), structural-homology
  (Foldseek 1-NN 0.876), classical ML, and other PLMs — i.e. the structure-aware embedding
  carries substrate-type information *beyond* homology (Stage 2).
- **The deployment model** is a 6-seed ProstT5 MLP (1024→512→256→4) with per-label
  MCC-optimal thresholds calibrated on **cluster-disjoint out-of-fold** predictions and an
  applicability-domain flag (Stage 2 → 3).
- **Focused beats broad.** A SAM-trained model out-performs a model trained on the broad
  201k-enzyme UP set on the SAM-MTase task (overall MCC 0.883 vs 0.824), driven mostly by
  the rare DNA class (Stage 2).
- **Tree-of-life sweep** of ~112k reference-proteome MTases (Stage 3) reveals
  taxonomic/lifestyle substrate trends (Stage 4) and a structured embedding landscape with
  recoverable, family-specific clusters (Stage 5).
- **Interpretation** (Stage 6): discriminative motifs recover the canonical catalytic /
  cofactor signatures, but the classifier's substrate signal localises to **buried fold-core
  residues**, not the conserved catalytic site — across all four substrate classes.

## Conventions used throughout

- **Embeddings.** ProstT5 (structure-aware, 1024-D) per-protein mean embeddings are the
  classifier features; ProtT5 embeddings are used **only for clustering** (homology splits),
  "to match the rest of the project". Embeddings are produced upstream by
  `Code/prostt5/embed_prostt5.py` into `Code/files/prostt5_*.h5` (not shipped here).
- **Labels.** Multi-label `[Other, Protein, RNA, DNA]`; single-substrate proteins are used
  for the headline analyses.
- **Validation.** Homology-aware **cluster-disjoint** cross-validation everywhere
  (greedy cosine clustering; SAM @0.6, UP @0.5). Per-label decision thresholds are
  **MCC-optimal, tuned on TRAIN/out-of-fold** predictions — never a fixed 0.5.
- **Compute.** Jobs run via SLURM `run_*.sh` drivers (conda env `classifier_env`, GPU).

## Reproduction note

Large regenerable inputs (per-protein/per-residue embedding `*.h5`, the raw corpus/UniProt
FASTAs, mmseqs temp dirs) are **not** shipped — they are produced by the upstream embedding
step and the stage scripts. Models, result tables, figures, and the Stage-6 structures are
included.
