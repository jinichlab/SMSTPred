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
