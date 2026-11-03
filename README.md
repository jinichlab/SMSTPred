# Predicting methyltransferase substrate specificity across the tree of life with structure-aware protein language models

Analysis code and results for the manuscript. The project builds a classifier that
predicts the **substrate type** of a SAM-dependent methyltransferase (MTase) —
**DNA / RNA / Protein / Metabolite (small-molecule)** — directly from a protein's sequence via
a **structure-aware protein language model (ProstT5)**, validates it against
homology and classical-ML baselines, deploys it across the tree of life, and interprets
*what* the model reads to make each call.

## The pipeline (stages)

| Stage | Folder | What it does |
|---|---|---|
| 0 | [`0_uniprot_pull/`](0_uniprot_pull) | **Data refresh** (optional, needs network): pulls the current Swiss-Prot SAM-methyltransferase set (InterPro IPR029063), the broad enzyme set (UP) and the 1,259 reference proteomes through the UniProt REST API, pins the UniProt release, and reports what changed against the archived data. The stages below run on the archived Zenodo data; Stage 0 produces the inputs for a re-run on newer UniProt. |
| 1 | [`1_clustering/`](1_clustering) | Homology-aware **cluster-disjoint** partitioning of the training sets so cross-validation never leaks homologues: the archived paper used greedy cosine clustering; the refresh uses **MMseqs2 at 30 % identity** and draws 16 splits of random clusters (65 % train / 15 % validation / 20 % test). |
| 2 | [`2_training/`](2_training) | The **classifier bake-off** (embedding MLPs trained with BCE and early stopping vs MMseqs2 / Foldseek / profile-HMM baselines, fixed 0.5 threshold; the archived paper's loss/epoch tuning and classical-ML arm are still there), the **SAM-vs-UP** training-set comparison, and the final **deployment model**. |
| 3 | [`3_predictions/`](3_predictions) | Scores new proteins with the deployment model (per-label score, call, ensemble disagreement, applicability-domain flag). |
| 4 | [`4_phylo_trends/`](4_phylo_trends) | Phylogenetic / taxonomic trends in the predicted substrate usage across the tree of life. |
| 5 | [`5_umap/`](5_umap) | Global ProstT5 UMAP of all tree-of-life MTases + substrate-coloured maps and cluster-enrichment figures. |
| 6 | [`6_dark_cluster_annotation/`](6_dark_cluster_annotation) | Cross-validates the 33 uncharacterized ("dark", zero-SwissProt-anchor) clusters surfaced in Stage 5 against two independent methods — **CLEAN** (EC-centroid nearest neighbor) and **Horizyn** (reaction-retrieval) — to flag genuine discovery leads vs clusters mis-swept in from other chemistry. |
| 7 | [`7_motif_analysis/`](7_motif_analysis) | **Interpretability**: discriminative sequence motifs (STREME, FIMO) and their structural mapping onto experimental enzyme–cofactor–substrate complexes — *what* distinguishes the substrate classes. |

## Two ways to use this repository

* **Run the classifier on your own sequences:** [`ready_classifier/`](ready_classifier) is a small, separate tool: a FASTA in, the substrate type of each protein out, with a step-by-step walkthrough of the input
  format and every output column. It always uses the newest trained model of this repository (`git pull` updates it).
* **Reproduce or refresh the whole analysis:** `bash workflow/run_refresh.sh` (`make refresh`) runs the entire paper as one SLURM job chain, from the UniProt pull to every figure and table, for the newest
  UniProt release; `bash workflow/setup_tools.sh` checks the external tools it needs ([`EXTERNAL_TOOLS.md`](EXTERNAL_TOOLS.md)). The outputs of a run are listed in `<release>/FIGURES.md` and the numbers
  quoted in the Results in `<release>/key_numbers.md`. `make verify` re-runs the analyses of the archived paper from the Zenodo data.

Each folder has its own `README.md` with scripts, inputs, outputs, and parameters. Composite-figure
and rendered-table builders are in [`figures/`](figures); the small data tables (training labels,
organism metadata) are in [`data/`](data).

## Refreshing the analysis with a new UniProt release

```bash
conda activate smstpred                 # make env; Foldseek: make tools; the embedders: make embeddings-env
make refresh-dry-run                    # prints the SLURM job chain and its dependencies, submits nothing
make refresh                            # pull -> bin -> cluster (+ SAM and UP splits) -> structures -> embeddings -> bake-off, classical ML, SAM-vs-UP (16-task arrays) -> reports
```

Each stage is its own job and starts only if the jobs before it finished without error. The settings that shape the results (seed, split fractions, epochs, early-stopping patience, threshold) are in
[`refresh_config.json`](refresh_config.json); the splits are drawn by the code from each release's clusters (same data and seed give the same splits). Not automatic: structures that AlphaFold DB cannot supply
need ColabFold (`embeddings/README.md`). Results land in `<release>/benchmark/` (Fig 2A, Table 1, Fig S4, paired tests, train/validation loss figures) and `<release>/sam_vs_up/` (Fig 2B). The numbers in "Headline results" below are those of the archived
Zenodo data (`make verify` reproduces them); the stages after the bake-off (3-7) are not part of the refresh chain yet. Details: [`workflow/run_refresh.sh`](workflow/run_refresh.sh), [`2_training/README.md`](2_training/README.md).

## Headline results

- **Embedding MLPs are statistically tied at the top.** Under cluster-disjoint CV (16 iterations) the ProstT5
  MLP reaches **MCC 0.926**, tied with the ProtT5+ESM-IF (0.928), SaProt (0.927), ProstT5+ESM-IF (0.920) and
  ProtT5 (0.919) MLPs; sequence homology (MMseqs2 1-NN 0.882) and structural homology (Foldseek 1-NN
  0.873) are not significantly separated either, while the ESM-IF-only MLP (0.859) and the profile-HMM (0.582,
  the floor) are significantly worse. On the same ProstT5 features a linear SVM (0.925) matches the MLP, so the
  representation, not the classifier head, drives performance. ProstT5 is deployed because it needs only
  sequence as input (Stage 2).
- **The deployment model** is a 6-seed ProstT5 MLP (1024→512→256→4) with per-label
  MCC-optimal thresholds calibrated on **cluster-disjoint out-of-fold** predictions and an
  applicability-domain flag (Stage 2 → 3).
- **Focused training is not worse than broad training.** A SAM-trained model scores overall MCC 0.908
  against 0.862 for one trained on the broad 201k-enzyme UP set, on held-out SAM-MTases; no individual
  substrate class differs significantly under the Nadeau-Bengio corrected paired test (the rare DNA class
  shows the largest gap, 0.913 vs 0.770) (Stage 2).
- **Tree-of-life sweep** of ~112k reference-proteome MTases (Stage 3) reveals
  taxonomic/lifestyle substrate trends (Stage 4) and a structured embedding landscape with
  recoverable, family-specific clusters (Stage 5), including 33 uncharacterized clusters
  with no reviewed anchor — two independent orthogonal methods flag several as genuine
  novel-methylation discovery leads and others as likely mis-assigned non-MTases (Stage 6).
- **Interpretation** (Stage 7): STREME recovers the canonical catalytic and cofactor
  signatures of each class (26 of 32 motifs class-specific; 22 of those 26 within 4 Å of the
  substrate or cofactor), yet about 12 % of enzymes carry no such motif and sit scattered
  within motif-bearing families; these enzymes are enriched for out-of-domain and fragment
  sequences, so a short discriminative motif is not what carries their calls, but this
  analysis does not show where their signal lies.

## Conventions used throughout

- **Embeddings.** ProstT5 (structure-aware, 1024-D) per-protein mean embeddings are the
  classifier features; ProtT5 embeddings are used **only for clustering** (homology splits),
  "to match the rest of the project". Embeddings are produced upstream by the model
  scripts (mean-pooled per protein); the precomputed embedding files are archived on Zenodo
  (see [`data/README.md`](data/README.md)), not shipped here.
- **Labels.** Multi-label `[Metabolite, Protein, RNA, DNA]`; single-substrate proteins are used
  for the headline analyses.
- **Validation.** Homology-aware **cluster-disjoint** cross-validation everywhere
  (greedy cosine clustering; SAM @0.6, UP @0.5). Per-label decision thresholds are
  **MCC-optimal, tuned on TRAIN/out-of-fold** predictions — never a fixed 0.5.
- **Compute.** GPU jobs have SLURM `run_*.sh` drivers (`sbatch run_x.sh`, or `bash run_x.sh` on any GPU machine); they use the active conda environment, else the `smstpred` environment that `make env` builds from `environment.yml`.

## Reproduction

```bash
git clone https://github.com/jinichlab/SMSTPred.git && cd SMSTPred
make env && conda activate smstpred                    # conda environment from environment.yml
make data DATA_SRC=zenodo:10.5281/zenodo.23655998      # the archived embeddings and tables, checksum-verified, into data/
make verify                                            # re-runs every CPU analysis in a scratch copy and compares with the committed results
make verify-gpu                                        # optional, on a GPU: re-trains the models and compares outcomes (stated tolerances)
```

`make help` lists the targets. [`REPRODUCING.md`](REPRODUCING.md) has the pass / fail table of the last run, the tolerances, and what is
not covered (external installs, multi-hour sweeps). The `run_*.sh` GPU drivers run under SLURM or as plain `bash`, in the environment above.

Large inputs (per-protein embedding `*.h5` files, the raw iTOL tree export) are archived on
**Zenodo** (doi:10.5281/zenodo.23655998; reserved, resolves once the record is published) rather than in git; per-residue embeddings, the raw
corpus/UniProt FASTAs and mmseqs temp dirs are regenerable and not shipped. Small data tables
are in [`data/`](data). Models, result tables, figures, and the Stage-6 structures are
included. Scripts resolve their inputs relative to the repository (data directory: `data/`, or
`$SMSTPRED_DATA`); [`data/README.md`](data/README.md) lists what is archived and what is not.

## License, citation and funding

- **Code:** MIT ([`LICENSE`](LICENSE)).
- **Data, result tables, figures and trained model weights** in this repository and in the Zenodo record:
  CC BY 4.0. The label tables and structures derive from UniProtKB (with Rhea and ChEBI annotation) and the
  AlphaFold Protein Structure Database, all CC BY 4.0, so please credit these sources as well when reusing them.
- **Citation:** see [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository"); the data record is
  doi:10.5281/zenodo.23655998.
- **Funding:** the authors thank the HHMI Hanna Gray Fellowship (GT16787) and the UC San Diego NIH FIRST Grant
  (NCI 1UA54CA272220) for funding.
