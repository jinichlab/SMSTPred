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
| `bakeoff.py`, `run_bakeoff.sh`, `test_bakeoff.py` | The **refreshed classifier bake-off**, one cluster-disjoint train / validation / test split per run (16-task SLURM array): six embedding MLPs (BCE, 2,000 epochs, early stopping keeps the best-validation weights) and the MMseqs2 / profile-HMM / Foldseek baselines, scored on the test clusters at a fixed 0.5 threshold. Reads only the release folder and `../refresh_config.json`. See "The refreshed bake-off" below. |
| `bakeoff_report.py`, `run_bakeoff_report.sh` | Merges the 16 splits; runs `classifier_mcc_stats.py`, `plot_cluster_mcc_zoom.py` (Fig 2A), `build_metrics_table.py` (Table 1) on that folder; the paired tests against the search baselines; the train / validation loss figures and the epochs early stopping kept. |
| `sam_vs_up.py`, `sam_vs_up_report.py`, `run_sam_vs_up.sh`, `run_sam_vs_up_report.sh`, `test_sam_vs_up.py` | **SAM-trained vs UP-trained ProtT5 MLP** on the same held-out SAM test clusters (Fig 2B of the refresh): 16 splits, the UP model trained on the UP clusters that remain after those containing a test protein are removed, early-stopped on 15 % of the UP clusters; paired tests incl. Nadeau-Bengio. See "SAM vs UP" below. |
| `archive_epoch_loss_selection/` | Retired tools that chose the epoch and the loss from a sweep (no longer used; why: its README). |
| `classifier_mcc_stats.py` | Friedman omnibus, Nemenyi / critical-difference diagram, Bonferroni-Dunn vs MLP_ProstT5 and Wilcoxon tests over the bake-off's cluster-disjoint MCCs → `stats/`. |
| `plot_cluster_mcc_zoom.py` | Zoomed cluster-disjoint MCC bar charts from the bake-off summary, bars coloured by method family: **Fig 2A** (the six embedding MLPs vs the three search-tool baselines), each with Friedman + Bonferroni-Dunn. |
| `build_metrics_table.py` | **Table 1**: every metric (MCC / AUPRC / AUROC / F1 / accuracy) for the nine representation and search-tool methods, rendered from the bake-off summary CSV. |
| `sam_vs_up_prostt5_comparison.py` | The **archived paper's** SAM-trained vs UP-trained ProstT5 MLP, both scored on the *same* held-out SAM test set (16 repeated cluster-disjoint 90/10 splits), per-label MCC-optimal thresholds. The SAM model is the deployed configuration (300 epochs); the UP model is trained for 400 epochs with the per-label losses that are best at that epoch in the UP sweep. Answers whether broad (UP) training helps the SAM-MTase task. `sam_vs_up_prostt5_stats.py` (Nadeau-Bengio corrected paired test) and `plot_sam_vs_up_prostt5.py` make Fig 2B; `plot_fig2_main.py` composes Fig 2. |
| `train_deployment_ensemble.py` | Builds the **deployment model**: frozen per-label thresholds from cluster-disjoint **out-of-fold** predictions + a 6-seed deep ensemble trained on all data + an applicability-domain cutoff. |
| `run_*.sh` | SLURM drivers. |

## The deployment model (`deployment_model/`)

ProstT5 MLP `1024→512→256→4`, **6-seed ensemble**, 300 epochs, per-label loss
`{Metabolite: sqrt_weighted_bce; Protein/RNA/DNA: bce}`. Files: `ensemble_member_{0..5}.pt`,
`thresholds.json` (OOF MCC-optimal cutoffs ≈ Metabolite 0.63 / Protein 0.63 / RNA 0.83 /
DNA 0.11), `model_config.json` (incl. applicability-domain cosine cutoff ≈ 0.94),
`oof_predictions.csv`. Consumed by Stage 3 `predict.py`.

## Settings behind the benchmark

- **Epochs and per-label losses.** Every full-set MLP trains for **300 epochs**; for each embedding and
  substrate label the loss (BCE, weighted BCE, sqrt-weighted BCE or focal) is the one with the highest mean
  MCC at epoch 300 in the embedding sweeps (Fig S2). An exact tie goes to the first loss in that order. The
  one exception is MLP_ProstT5, which keeps the deployed loss map (Metabolite: sqrt-weighted BCE) so that the
  benchmark scores the model that ships; the rule alone would pick weighted BCE there, a 0.001-MCC tie. The
  UP-trained ProstT5 in Fig 2B uses 400 epochs and the per-label best losses at epoch 400 in the UP sweep (Fig S5).
  All of this is in the bake-off settings and `UP_CFG`
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
  Δ = +0.001 MCC, SaProt better in 8 of 16 iterations, Wilcoxon p = 0.90; Holm-adjusted p = 1.00 across the
  eight comparisons of `paired_tests_vs_control.csv`, written by `benchmark_supplement_stats.py`).
- **Per-label loss** matters mainly for the rare DNA class; abundant classes are
  loss-insensitive (the MCC-curve grids show all losses converging).
- **SAM-vs-UP**: focused SAM training is not worse than broad UP training on held-out SAM (overall MCC
  0.908 vs 0.862); no individual substrate class differs significantly under the Nadeau-Bengio corrected
  paired test (two-sided p = 0.65 Metabolite, 0.93 Protein, 0.75 RNA, 0.25 DNA, 0.30 overall), the largest gap
  being **DNA** (0.913 vs 0.770), which a few splits drive. The UP model even had a mild homology advantage
  (only exact test entries removed).

## Outputs

`mlp_embedding_comparison_output/` (+ `mlp_embedding_comparison_saprot_output/` for SaProt),
`sam_classifier_comparison_output/` (bake-off CSVs for the full set, `stats/`, Fig 2A, Fig S4, Fig 2, Table 1;
`PROVENANCE.md` says how it was assembled), `UP_prostt5_loss_comparison_output/`,
`sam_vs_up_prostt5_comparison_output/` (Fig 2B), `train_val_loss_output/` (Figs S3, S6, summary CSVs),
`deployment_model/` — summary/detailed/per-label CSVs, MCC grids, the zoom bar charts, the SAM-vs-UP figure,
and the trained model.

External tools used by the bake-off: `mmseqs`, `foldseek`, `mafft`, `hmmbuild`/`hmmsearch`.


## The refreshed bake-off (`bakeoff.py`)

The protocol, fixed in [`../refresh_config.json`](../refresh_config.json) and written into every output folder (`config_used.json`):

1. **Splits** (Stage 1, `1_clustering/tvt_splits.py`): 16 times, 20 % of the 30 %-identity MMseqs2 clusters are the **test** set, 15 % the **validation** set, 65 % the **train** set. A method never sees the clusters of
   the other two parts.
2. **Embedding MLPs** (ProstT5, ProtT5, SaProt, ESM-IF, ProtT5+ESM-IF, ProstT5+ESM-IF; in -> 512 -> 256 -> 4): **plain BCE** on all four labels, Adam (lr 0.001), one step per epoch on a random 90 % of the train proteins,
   trained for **all 2,000 epochs** so the train / validation loss can be plotted to the end. **Early stopping** (patience 80 epochs) only decides which weights are kept: those of the epoch with the lowest total
   validation loss before the loss failed to improve for 80 epochs in a row. The test proteins are predicted with those weights.
3. **Decision threshold: a fixed 0.5** for the MLPs and for the MMseqs2 / Foldseek vote (a bit-score-weighted majority vote over the top-10 hits' labels, searched against the train proteins only). The profile HMM scores are bit
   scores, which have no 0.5: its per-label threshold is the MCC-optimal one on the validation scores. `mcc_valtuned` records, for every method, the MCC at per-label thresholds tuned on the validation predictions: on the 2026_03
   release that was *lower* than 0.5 for every method (a validation set holds 9-33 DNA and at least 40 Metabolite / 61 Protein positives, too few for a stable optimum), so 0.5 is the primary choice.
4. **Statistics**: Friedman / Nemenyi / Bonferroni-Dunn / Wilcoxon over the 16 splits (`classifier_mcc_stats.py`) and `paired_tests_vs_baselines.csv` (mean difference, bootstrap interval, splits won, Wilcoxon p, Holm over
   each MLP's three comparisons). The 16 splits are random re-draws of the same proteins and overlap, so the p-values are optimistic: read the differences and the split counts first.

Why this design (the choices made on the 2026_03 data, in short): choosing the epoch or the loss from a sweep scored on the held-out clusters is not independent of the test set; the per-label loss choice was within noise
(`archive_epoch_loss_selection/`); thresholds tuned on the *training* predictions were unstable (a trained network is near-perfect on its training data) and cost 0.01-0.04 MCC against 0.5; plain BCE was as good as any weighted
loss, also for DNA, whose ranking is already near-perfect (AUROC 0.9985).

```bash
make splits          # (also part of `make cluster`)   make bakeoff   # sbatch array, ~10 min per split   make bakeoff-report   # merge, statistics, figures
make bakeoff-test    # offline tests: early stopping, the 0.5 / validation-tuned thresholds, the loader, the splits, the threshold search
```
Release 2026_03, 16 splits, mean test MCC (fixed 0.5; HMM validation-tuned): ProtT5 0.940, SaProt 0.938, ProtT5+ESM-IF 0.937, ProstT5 0.895, Foldseek 0.887, ProstT5+ESM-IF 0.883, MMseqs2 0.848,
ESM-IF 0.823, profile HMM 0.519 (Friedman p = 3e-15). Early stopping kept the epoch (median over the splits) at 204 for ProtT5, 215 ProstT5, 236 ProtT5+ESM-IF, 258 ProstT5+ESM-IF, 308 ESM-IF and 336 SaProt, with
a wide spread between splits (e.g. ProstT5 90-486); it fired in every split. Re-running the same code reproduced every MLP, MMseqs2 and Foldseek number exactly; the profile HMM varied by up to 0.045 MCC in a split
between two identical runs (cause not found; MAFFT was single-threaded in both).

Output folder `<release>/benchmark/`: `iterations/iterNN_*` (per split: metrics, per-label metrics, the train / validation loss of every epoch, the epoch kept, the validation and test scores of every method), the merged
`sam_classifier_comparison_*.csv`, `stats/`, `sam_classifier_cluster_mcc_zoom.png` (Fig 2A), `metrics_table.png` (Table 1), `paired_tests_vs_baselines.csv`, `train_val_loss_by_embedding.png`,
`train_val_loss_by_label.png`, `kept_epochs.csv`, `config_used.json`.

### SAM vs UP (`sam_vs_up.py`)

Does a model trained on the broad UP enzyme set (about 200,000 proteins) predict the substrate type of held-out SAM methyltransferases as well as one trained on SAM? Per split `i`: the **test** set is the 20 % of
SAM clusters held out in the bake-off; every UP cluster (UP and SAM clustered together at 30 %, `mmseqs_clusters.py up`) that contains one of those proteins is **removed** from the UP set
(`1_clustering/up_splits.py`); 15 % of the remaining UP clusters are the **validation** set, the rest the **train** set (about 172,600 / 29,700 proteins). Both models are ProtT5 MLPs trained the same way (BCE,
2,000 epochs, best validation epoch kept, patience 80), both scored at a fixed 0.5 on the same test proteins. Arms: `ProtT5_SAM` (the bake-off's `MLP_ProtT5`, not retrained), `ProtT5_UP`, and a sensitivity arm
`ProtT5_SAMclean`: the SAM MLP retrained without the SAM train / validation proteins that share a joint SAM+UP cluster with a test protein (the UP arm loses them, the SAM arm of the bake-off keeps them: 511-1,518
proteins per split, because SAM's own clusters are made on SAM alone). The tests are two-sided: Wilcoxon, the paired t-test and the **Nadeau-Bengio corrected resampled t-test** (variance inflated by 1/n + n_test/n_train,
here 0.334; the splits are cluster draws, which the formula does not model, so it is an approximation); the first two ignore that the training sets of the 16 splits overlap and are optimistic.
Release 2026_03, MCC (mean over 16 splits): SAM 0.940, UP 0.911, SAMclean 0.918; UP minus SAM -0.030 (SAM better in 14/16 splits; Wilcoxon p = 0.016, paired t 0.038, Nadeau-Bengio 0.38); by label the gap is
DNA -0.065 (UP never better, 0/16; Wilcoxon p = 0.001, Nadeau-Bengio 0.15), Protein -0.046, Metabolite -0.002, RNA -0.005. Against the SAMclean arm the overall gap is -0.007 (not significant); the DNA gap stays
(UP 0.898, SAMclean 0.966). SAMclean also has fewer training proteins (about 4,100 against 4,900), so its drop from SAM mixes less homology with less data.

### The deployment model of the refresh (`deploy_ensemble.py`)

A **6-member ProtT5 MLP ensemble**; each member is the bake-off's MLP (BCE, 2,000 epochs trained, best validation epoch kept, patience 80) trained on a different 5/6 of the clusters of the deployment set
(`in_deployment`, <= 2,000 aa, 7,582 proteins in 492 clusters) and early-stopped on the other 1/6: a 6-fold **cluster partition** packed to balance protein count and the positives of every label. The prediction is
the mean probability at a fixed 0.5 per label; the spread across members is the disagreement signal. Why this and not the archived "6 seeds on all the data, 300 epochs": on the 16 held-out SAM test sets (nested: members
trained on the 80 % non-test clusters) six seeds on the same data gave 0.9409 against 0.9405 for one model (nothing), six members with different cluster splits 0.9558 (+0.015, better in 15 of 16 test sets, all four labels;
Nadeau-Bengio p 0.07-0.09), and a partition was as accurate as six independent random 15 % validation draws (0.9558 both) while validating every protein once (random draws leave 38 % never validated); the gain
saturates at about four members.

Outputs in `<release>/deployment/` are the files `3_predictions/predict.py` reads, so it is a drop-in (`DEPLOY_DIR=<release>/deployment`): `ensemble_member_{0..5}.pt`, `model_config.json`, `thresholds.json` (0.5),
`ad_reference.h5` (the training vectors) with the **applicability-domain cutoff** (1st percentile of the training proteins' nearest-neighbour cosine similarity: 0.813 in ProtT5 space), plus `oof_predictions.csv`,
`calibration.csv`, `kept_epochs.csv`, `curves.csv.gz` and the loss figure. **Out-of-fold predictions come from separate audit members**: a deployed member chose its stopping epoch on its own validation fold, so its
predictions there are optimistic (release 2026_03: MCC 0.974, against 0.955 for audit members whose test fold was never used for training or for the epoch; the optimism does not depend on the number of folds).
The honest figure, a single member trained on 67 % of the clusters, is consistent with the nested test (0.946-0.949 per member on unseen clusters; the ensemble 0.956).

**ProtT5 vectors of the proteomes.** The training vectors (half-precision encoder, mean over the residue tokens) and the archived tree-of-life vectors (full-precision model, mean over every non-padding token,
end-of-sequence token included) come from different pipelines. Checked on the 1,104 proteins in both: median cosine 1.0000, median relative L2 difference 0.7 % (0.4 % at 500-1,000 residues: the end-of-sequence token), a systematic
shift of 0.6 % of a vector norm; the deployment ensemble gives the **same calls for all 1,104 proteins** (every label) and the same applicability-domain flag, mean |probability difference| 0.0002 (max 0.11, one protein
whose vector differs, cosine 0.88; two of the four proteins with a cosine below 0.999 had their sequence modified after the tree-of-life files were made, for the other two I cannot tell which UniProt version those files used).
