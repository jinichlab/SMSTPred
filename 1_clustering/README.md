# Stage 1 — Homology-aware clustering & threshold selection

Builds the **cluster-disjoint** partitions used by every downstream cross-validation, so
that a held-out protein never shares a near-homologue with the training pool. Clustering is
**greedy single-linkage on cosine similarity** of the ProtT5 embeddings (ProtT5 is used for
clustering throughout the project, even though ProstT5 is the classifier feature).

## Scripts

| Script | Purpose |
|---|---|
| `sam_clustering.py` | Cluster the SAM_swiss MTase set and justify the **0.6** cosine threshold. Sweeps thresholds, computes within-set vs between-set similarity baselines, and writes the cluster assignment used by Stages 2–3. |
| `up_clustering.py` | Same for the broad **UP_swiss** enzyme set (201,339 proteins); the looser **0.5** threshold is justified by the within-set similarity distribution. Drives the SAM-vs-UP comparison. |
| `run_clustering.sh`, `run_sam_clustering.sh` | SLURM drivers. |

## Key parameters / numbers

- **SAM clustering: threshold 0.6** → ~45 clusters over 7,405 proteins (one large
  cluster dominates; clusters are packed into folds balanced by protein count).
- **UP clustering: threshold 0.5** → ~406 clusters over 201,339 proteins (loose threshold,
  large clusters; full pairwise clustering at scale is avoided — see the script notes).

## Outputs

- `sam_clustering_output/`, `up_clustering_output/` — cluster pickles (`*_clusters_t*.pkl`,
  reused downstream so the partition is identical across stages), threshold-selection
  figures and similarity-baseline CSVs.

## Why it matters

Substrate type tracks protein family, so a *random* split would let the classifier
memorise families and over-report accuracy. Cluster-disjoint splits force genuine
generalisation to **remote homologues**, which is the regime deployment actually faces.

## MMseqs2 clustering of the typed sets (the refresh, runs right after the substrate typing)

`mmseqs_clusters.py` (driver `run_mmseqs_clusters.sh`, `make cluster`) clusters **by sequence** instead of by embedding, so it needs no embeddings and the
cluster-disjoint splits are fixed as soon as the labelled set is known, before anything is embedded. It reads `labels_SAM.csv` / `labels_UP.csv` of the newest
pull (substrate typing: `0_uniprot_pull/bin_substrates.py`) and the sequences of the pulled tables.

| Run | Proteins | Used for |
|---|---|---|
| `sam --subset benchmark` | the SAM proteins with `in_benchmark` (<= 1,000 aa) | Fig 2A, Table 1, S2 / S3, every cross-validation |
| `sam --subset deployment` | `in_deployment` (<= 2,000 aa) | the deployment model, its out-of-fold thresholds |
| `up` | UP `in_benchmark` + the SAM benchmark proteins (clustered together, so UP homologues of a held-out SAM protein can be removed from UP training) | Fig 2B |

MMseqs2 easy-cluster at **30 % identity**, 80 % of both sequences aligned (`--cov-mode 0`), `-s 7.5`, greedy set cover, `--cluster-reassign 1`: the settings of the side project
`~/Code/mmseqs_clustering`. Checked against it: on the identical 7,405 archived proteins the port gives the identical partition (adjusted Rand index 1.000000, the same cluster id
for every protein) and 16 identical train / validation splits. Splits: shuffle the clusters, first 90 % of the clusters train, the rest validation, redraw until every label is in
validation; 16 draws from the seed-42 stream (the first 6 are the S2 / S3 splits). There is no fold column: the benchmark uses these repeated cluster-level draws, not K folds.

Release 2026_03: 7,553 benchmark proteins in 502 clusters (248 singletons; the ten largest hold 43.9 % of the proteins, the largest 830 RNA enzymes), 7,582 deployment proteins in 492.
Outputs go to `<release>/clusters/` (`clusters_<tag>.csv`, `splits_<tag>.json`, summaries, raw MMseqs2 files and parameters) and two columns, `cluster_benchmark` /
`cluster_size_benchmark` (or `_deployment`), are written into the label table; re-running `bin_substrates.py` rewrites the table without them, so cluster after typing.
The ProtT5-cosine scripts above stay for reproducing the archived paper (their pickles are for the archived 7,405 proteins).

## The train / validation / test splits of the bake-off (`tvt_splits.py`, runs right after the SAM benchmark clustering)

`tvt_splits.py` (called by `run_mmseqs_clusters.sh`, or `make splits`) turns the `cluster_benchmark` column into the **16 splits of the classifier bake-off**: for every split the clusters are
shuffled and **20 % of the clusters are the test set, 15 % the validation set and the remaining 65 % the training set** (`round(0.20 * 502) = 100`, `75`, `327` clusters in release 2026_03).
Settings (fractions, number of splits, seed 42) are in [`../refresh_config.json`](../refresh_config.json); the same data and seed always give the same splits, a new release (new clusters) gives new ones.
Nothing is balanced or redrawn: what a partition holds is what chance gives. Cluster sizes are very uneven (248 of the 502 clusters are singletons, the largest holds 830 proteins), so 20 % of the
clusters is 16-32 % of the proteins; the table printed by the script (and `splits_tvt_summary_*.csv`) lists proteins and positives of every label per split and partition. Release 2026_03: every test set
has at least 15 DNA, 53 Metabolite, 70 Protein and 640 RNA positives; the validation sets at least 9 DNA (one split), 40, 61 and 200.
Output: `<release>/clusters/splits_tvt_sam_benchmark_id30_t20_v15.json` (`train` / `validation` / `test` lists of UniProt entries per iteration; the fractions are part of the name) and a summary CSV.
The older 90 / 10 train / validation splits written by `mmseqs_clusters.py` stay (the SAM-vs-UP comparison and the archived figures use them); the bake-off does not.

## The UP train / validation sets of the SAM-vs-UP comparison (`up_splits.py`, `make up-splits`)

For every SAM test set of the bake-off, `up_splits.py` removes from the UP set **every UP cluster that contains a test protein** (the UP clustering, `mmseqs_clusters.py up`, clusters UP and the SAM benchmark
proteins together at 30 %, so this takes out the test proteins and their UP homologues), then draws **15 %** of the remaining UP clusters as the validation set (seeded by the split number) and keeps the rest as
the train set. Release 2026_03 per split: 168,000-175,000 train and 26,000-34,000 validation proteins, at least 175 DNA / 2,700 Protein positives in a validation set; 105-123 clusters are removed with the test
proteins. UP-only proteins are removed in only 6 of the 16 splits (1-30 each, 36 in all); the bulk of what goes with the test clusters (511-1,518 proteins per split) are SAM train / validation proteins that
share a joint cluster with a test protein. Output: `<release>/clusters/splits_tvt_up_id30_t20_v15.json` and a summary CSV. Needs `CLUSTER_WHAT=up` clustering first (about 210,000 sequences: 32 cpus, 100 GB).
