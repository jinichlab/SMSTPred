# Stage 2 — Model selection, benchmarking & the deployment model

Everything between "we have embeddings + cluster-disjoint splits" and "we have a deployed
classifier": the per-label loss/epoch tuning, the head-to-head classifier benchmark, the
training-set (SAM-vs-UP) comparison, and the final production model.

## Scripts

| Script | Purpose |
|---|---|
| `mlp_loss_comparison.py` | Per-label **loss-function & training-epoch** tuning under cluster-disjoint CV, for SAM / SAM-balanced / UP (ProtT5). Produces the per-label "which loss wins" CSVs that motivate the per-label mixed loss. |
| `mlp_loss_comparison_prostt5.py` | The **ProstT5** version of the loss comparison for the UP set (so SAM-vs-UP is apples-to-apples on the deployment embedding). Dense-early epoch grid (ProstT5 peaks fast), 6 cluster-disjoint iterations. |
| `mlp_embedding_comparison.py`, `mlp_embedding_comparison_saprot.py` | The **embedding × loss × epoch sweeps** behind the MCC-vs-epoch grids for the SAM full and balanced sets (Figs S2/S3): ProtT5, ProstT5, ESM-IF and the ProtT5+ESM-IF / ProstT5+ESM-IF concatenations (`..._saprot.py`: SaProt), 4 loss types, 6 cluster-disjoint splits, 20,000 epochs with dense-early checkpoints. Write the per-label CSVs that `plot_mcc_curves_grid.py` reads. |
| `plot_mcc_curves_grid.py` | Per-label **MCC-vs-epoch grids** (one figure family): SAM full, SAM balanced (rows = embedding), and UP (single ProstT5 row); lines = loss types, ±1 SD bands. |
| `sam_full_vs_balanced_classifier_comparison.py` | The **classifier bake-off**: 16 methods — MLPs over ProtT5 / ProstT5 / ESM-IF / SaProt / concatenations, homology baselines (**MMseqs2 1-NN**, **Foldseek 1-NN**, **Profile-HMM**), and classical ML (XGBoost, RF, SVMs, LogReg) — on a **shared held-out test set** (cluster + random), full vs class-balanced. Per-class MCC-optimal thresholds for every method. |
| `plot_cluster_mcc_zoom.py` | Zoomed cluster-disjoint MCC bar chart from the bake-off summary (bars coloured by method family). |
| `build_metrics_table.py` | **Table 1**: every metric (MCC / AUPRC / AUROC / F1 / accuracy) for the nine representation and search-tool methods, rendered from the bake-off summary CSV. |
| `sam_vs_up_prostt5_comparison.py` | **SAM-trained vs UP-trained** ProstT5 MLP, both scored on the *same* held-out SAM test set (16 repeated cluster-disjoint 90/10 splits), per-label MCC-optimal thresholds. Answers whether broad (UP) training helps the SAM-MTase task. |
| `train_deployment_ensemble.py` | Builds the **deployment model**: frozen per-label thresholds from cluster-disjoint **out-of-fold** predictions + a 6-seed deep ensemble trained on all data + an applicability-domain cutoff. |
| `run_*.sh` | SLURM drivers. |

## The deployment model (`deployment_model/`)

ProstT5 MLP `1024→512→256→4`, **6-seed ensemble**, 300 epochs, per-label loss
`{Other: sqrt_weighted_bce; Protein/RNA/DNA: bce}`. Files: `ensemble_member_{0..5}.pt`,
`thresholds.json` (OOF MCC-optimal cutoffs ≈ Other 0.63 / Protein 0.63 / RNA 0.83 /
DNA 0.11), `model_config.json` (incl. applicability-domain cosine cutoff ≈ 0.94),
`oof_predictions.csv`. Consumed by Stage 3 `predict.py`.

## Key results

- **MLP_ProstT5 is statistically tied for best** (cluster-disjoint full-set MCC **0.928**;
  MLP_SaProt 0.934), with MMseqs2 1-NN 0.882 (sequence homology, not significantly
  separated), ProtT5 0.875 and Foldseek 1-NN 0.873 (structural homology) just behind and
  the Profile-HMM the floor (0.582). A linear SVM on the same ProstT5 features (0.925)
  matches the MLP, so the representation, not the classifier head, drives performance.
  ProstT5 is deployed because it needs only sequence as input. All means over 16
  cluster-disjoint iterations.
- **Per-label loss** matters mainly for the rare DNA class; abundant classes are
  loss-insensitive (the MCC-curve grids show all losses converging).
- **SAM-vs-UP**: focused SAM training is at least as good as broad UP training on held-out
  SAM (overall MCC 0.908 vs 0.826); no individual substrate class differs significantly
  under the Nadeau-Bengio corrected paired test, the largest gap being **DNA** (0.913 vs
  0.729). The UP model even had a mild homology advantage (only exact test entries removed).

## Outputs

`mlp_embedding_comparison_output/` (+ `mlp_embedding_comparison_saprot_output/` for SaProt), `sam_classifier_comparison_output/`,
`UP_prostt5_loss_comparison_output/`, `sam_vs_up_prostt5_comparison_output/`,
`deployment_model/` — summary/detailed/per-label CSVs, MCC grids, the zoom bar chart, the
SAM-vs-UP figure, and the trained model.

External tools used by the bake-off: `mmseqs`, `foldseek`, `mafft`, `hmmbuild`/`hmmsearch`.
