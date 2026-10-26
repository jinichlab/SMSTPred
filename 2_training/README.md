# Stage 2 — Model selection, benchmarking & the deployment model

Everything between "we have embeddings + cluster-disjoint splits" and "we have a deployed
classifier": the per-label loss/epoch tuning, the head-to-head classifier benchmark, the
training-set (SAM-vs-UP) comparison, and the final production model.

## Scripts

| Script | Purpose |
|---|---|
| `mlp_loss_comparison.py` | Per-label **loss-function & training-epoch** tuning under cluster-disjoint CV, for SAM / UP (ProtT5; the script also covers the balanced SAM subset, which the manuscript does not use). Produces the per-label "which loss wins" CSVs that motivate the per-label mixed loss. |
| `mlp_loss_comparison_prostt5.py` | The **ProstT5** version of the loss comparison for the UP set (so SAM-vs-UP is apples-to-apples on the deployment embedding). Dense-early epoch grid (ProstT5 peaks fast), 6 cluster-disjoint iterations. |
| `mlp_embedding_comparison.py`, `mlp_embedding_comparison_saprot.py` | The **embedding × loss × epoch sweeps** behind the MCC-vs-epoch grids for the SAM full set (Fig S2): ProtT5, ProstT5, ESM-IF and the ProtT5+ESM-IF / ProstT5+ESM-IF concatenations (`..._saprot.py`: SaProt), 4 loss types, 6 cluster-disjoint splits, 20,000 epochs with dense-early checkpoints (the balanced subset is swept too; not used in the manuscript). Write the per-label CSVs that `plot_mcc_curves_grid.py` reads and from which the benchmark's epochs and per-label losses are chosen. |
| `plot_mcc_curves_grid.py` | Per-label **MCC-vs-epoch grids** (Figs S2, S5): SAM full (rows = embedding) and UP (single ProstT5 row); lines = loss types, ±1 SD bands. `--balanced` adds the balanced-subset grid. |
| `train_val_loss_curves.py`, `plot_train_val_loss_curves.py`, `run_train_val_loss.sh` | **Train vs validation loss per epoch** (Figs S3, S6): re-runs the same sweeps (same splits, MLP and 20,000-epoch budget) while recording the loss on the training pool and on the cluster-disjoint validation set, to show where validation loss bottoms out. The raw curves (about 50 MB) are on Zenodo (`train_val_loss_curves.tar`); extract them into `train_val_loss_output/`, then run the plotting script. SLURM array tasks 0–5 are SAM full, 12–15 UP; 6–11 (balanced subset) are optional. |
| `sam_classifier_comparison.py` | The **classifier bake-off**: 16 methods — MLPs over ProtT5 / ProstT5 / ESM-IF / SaProt / concatenations, homology baselines (**MMseqs2 1-NN**, **Foldseek 1-NN**, **Profile-HMM**), and classical ML (XGBoost, RF, SVMs, LogReg) — trained on the full SAM set and scored on a **shared held-out test set** (cluster + random splits, 16 iterations). Per-class MCC-optimal thresholds for every method. MLP epochs and per-label losses are set in `MLP_VARIANTS`. Environment options: `BAKEOFF_ONLY=mlp` (retrain only the MLPs; the splits stay identical), `BAKEOFF_MLP_VARIANTS=a,b`, `BAKEOFF_OUT`, and `BAKEOFF_DATASETS=full,balanced` (adds a class-balanced training arm; off by default, not in the manuscript). |
| `classifier_mcc_stats.py` | Friedman omnibus, Nemenyi / critical-difference diagram, Bonferroni-Dunn vs MLP_ProstT5 and Wilcoxon tests over the bake-off's cluster-disjoint MCCs → `stats/`. |
| `merge_bakeoff_rows.py` | Swaps re-trained MLP rows into the bake-off CSVs (how the committed results were assembled; see `sam_classifier_comparison_output/PROVENANCE.md`). |
| `plot_cluster_mcc_zoom.py` | Zoomed cluster-disjoint MCC bar charts from the bake-off summary, bars coloured by method family: **Fig 2A** (the six embedding MLPs vs the three search-tool baselines) and **Fig S4** (MLP_ProstT5 vs seven classical classifiers on the same ProstT5 features), each with Friedman + Bonferroni-Dunn. |
| `build_metrics_table.py` | **Table 1**: every metric (MCC / AUPRC / AUROC / F1 / accuracy) for the nine representation and search-tool methods, rendered from the bake-off summary CSV. |
| `sam_vs_up_prostt5_comparison.py` | **SAM-trained vs UP-trained** ProstT5 MLP, both scored on the *same* held-out SAM test set (16 repeated cluster-disjoint 90/10 splits), per-label MCC-optimal thresholds. The SAM model is the deployed configuration (300 epochs); the UP model is trained for 400 epochs with the per-label losses that are best at that epoch in the UP sweep. Answers whether broad (UP) training helps the SAM-MTase task. `sam_vs_up_prostt5_stats.py` (Nadeau-Bengio corrected paired test) and `plot_sam_vs_up_prostt5.py` make Fig 2B; `plot_fig2_main.py` composes Fig 2. |
| `train_deployment_ensemble.py` | Builds the **deployment model**: frozen per-label thresholds from cluster-disjoint **out-of-fold** predictions + a 6-seed deep ensemble trained on all data + an applicability-domain cutoff. |
| `run_*.sh` | SLURM drivers. |

## The deployment model (`deployment_model/`)

ProstT5 MLP `1024→512→256→4`, **6-seed ensemble**, 300 epochs, per-label loss
`{Other: sqrt_weighted_bce; Protein/RNA/DNA: bce}`. Files: `ensemble_member_{0..5}.pt`,
`thresholds.json` (OOF MCC-optimal cutoffs ≈ Other 0.63 / Protein 0.63 / RNA 0.83 /
DNA 0.11), `model_config.json` (incl. applicability-domain cosine cutoff ≈ 0.94),
`oof_predictions.csv`. Consumed by Stage 3 `predict.py`.

## Settings behind the benchmark

- **Epochs and per-label losses.** Every full-set MLP trains for **300 epochs**; for each embedding and
  substrate label the loss (BCE, weighted BCE, sqrt-weighted BCE or focal) is the one with the highest mean
  MCC at epoch 300 in the embedding sweeps (Fig S2). An exact tie goes to the first loss in that order. The
  one exception is MLP_ProstT5, which keeps the deployed loss map (Other: sqrt-weighted BCE) so that the
  benchmark scores the model that ships; the rule alone would pick weighted BCE there, a 0.001-MCC tie. The
  UP-trained ProstT5 in Fig 2B uses 400 epochs and the per-label best losses at epoch 400 in the UP sweep (Fig S5).
  All of this is in `MLP_VARIANTS` (`sam_classifier_comparison.py`) and `UP_CFG`
  (`sam_vs_up_prostt5_comparison.py`).
- **Train vs validation loss (Figs S3, S6).** Validation loss is lowest early (mean-curve minima at epochs
  35–357 on the full set, 104–411 for UP) and rises while training loss keeps falling, but MCC stays on its
  plateau, so epochs were chosen on MCC, with validation loss as a check: at 300 epochs the full-set
  validation loss is a median 1.2× its minimum, and at 400 epochs the UP loss is a median 1.1× its minimum.
  Spikes in the training loss (optimizer instability) start only after epoch 9,300 on the full set (weighted
  losses) and between epochs 2,260 and 5,470 for UP, all well beyond the chosen epochs.

## Key results

Cluster-disjoint CV, 16 iterations, full SAM training set; means over the iterations.

- **Five of the six embedding MLPs are statistically tied for best** (Friedman χ²(8) = 54.3, p = 6.0e-9
  across the nine representation and search-tool methods; Bonferroni-Dunn vs MLP_ProstT5): MCC
  ProtT5+ESM-IF **0.928**, SaProt 0.927, ProstT5 **0.926 (deployed)**, ProstT5+ESM-IF 0.920, ProtT5 0.919.
  Sequence homology (MMseqs2 1-NN 0.882) and structural homology (Foldseek 1-NN 0.873) are not significantly
  separated from them either; the ESM-IF-only MLP (0.859) and the profile-HMM (0.582, the floor) are
  significantly worse. ProstT5 is deployed because it needs only sequence as input (SaProt − ProstT5:
  Δ = +0.001 MCC, SaProt better in 8 of 16 iterations, Wilcoxon p = 0.90).
- **The representation, not the classifier head, drives performance.** On the same ProstT5 features the linear
  (0.925), RBF (0.921) and polynomial (0.916) SVMs are tied with the MLP (0.926); logistic regression,
  XGBoost and random forest are significantly worse, and the sigmoid-kernel SVM (0.876) is borderline (its
  rank gap sits just under the critical difference; Friedman χ²(7) = 65.5, p = 1.2e-11).
- **Per-label loss** matters mainly for the rare DNA class; abundant classes are
  loss-insensitive (the MCC-curve grids show all losses converging).
- **SAM-vs-UP**: focused SAM training is not worse than broad UP training on held-out SAM (overall MCC
  0.908 vs 0.862); no individual substrate class differs significantly under the Nadeau-Bengio corrected
  paired test (two-sided p = 0.65 Other, 0.93 Protein, 0.75 RNA, 0.25 DNA, 0.30 overall), the largest gap
  being **DNA** (0.913 vs 0.770), which a few splits drive. The UP model even had a mild homology advantage
  (only exact test entries removed). An earlier UP run with 3,500 epochs scored lower (overall 0.826) and
  the overall difference was then significant (p = 0.042); see the `UP_EPOCHS` note in the script.

## Outputs

`mlp_embedding_comparison_output/` (+ `mlp_embedding_comparison_saprot_output/` for SaProt),
`sam_classifier_comparison_output/` (bake-off CSVs for the full set, `stats/`, Fig 2A, Fig S4, Fig 2, Table 1;
`PROVENANCE.md` says how it was assembled), `UP_prostt5_loss_comparison_output/`,
`sam_vs_up_prostt5_comparison_output/` (Fig 2B), `train_val_loss_output/` (Figs S3, S6, summary CSVs),
`deployment_model/` — summary/detailed/per-label CSVs, MCC grids, the zoom bar charts, the SAM-vs-UP figure,
and the trained model.

External tools used by the bake-off: `mmseqs`, `foldseek`, `mafft`, `hmmbuild`/`hmmsearch`.
