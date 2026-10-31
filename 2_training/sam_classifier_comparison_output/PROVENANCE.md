# How the committed bake-off results were assembled

`sam_classifier_comparison_{detailed,per_label,summary}.csv` hold the classifier comparison on the **full SAM
training set** (cluster-disjoint and random splits, 16 iterations, one shared held-out test set per iteration).
Every figure, table and statistic in this folder (Fig 2A, Fig S4, Table 1, `stats/`) is computed from them.

## Where the rows come from

| Rows | Source |
|---|---|
| XGBoost, Random Forest, Logistic Regression, the four SVMs, MMseqs2 1-NN, Foldseek 1-NN, Profile-HMM | The original bake-off run: all methods trained from scratch with `sam_classifier_comparison.py` (then named `sam_full_vs_balanced_classifier_comparison.py`, commit `72c3d17`). |
| MLP over ProtT5, ESM-IF, ProtT5+ESM-IF, ProstT5+ESM-IF, SaProt | MLP-only re-run (`BAKEOFF_ONLY=mlp`, SLURM job 300280) with 300 epochs and, per embedding and label, the loss that is best at epoch 300 in the Fig S2 sweep: the settings in `MLP_VARIANTS`. Raw rows: `provenance/mlp_rerun_all_six_variants_full.csv`. |
| MLP over ProstT5 | MLP-only re-run of that variant alone with the **deployed** loss map (Metabolite: sqrt-weighted BCE; `BAKEOFF_MLP_VARIANTS=MLP_ProstT5`, SLURM job 300283). The first re-run had used weighted BCE for the Metabolite head, which the selection rule picks by a 0.001-MCC tie; the deployed map was kept so the benchmark scores the shipped model. Raw rows: `provenance/mlp_rerun_prostt5_deployed_full.csv`. |

The MLP settings of the original run (100–400 epochs and other loss maps) are in git history at `72c3d17`.
The original run also trained on a class-balanced subset; those rows were dropped, because the manuscript does not use
that arm (`BAKEOFF_DATASETS=full,balanced` still runs it). The rows were combined with `merge_bakeoff_rows.py`
(`--keep-datasets full`).

## Why combining runs is valid

The train/test split of every iteration is a function of the random stream, and the Profile-HMM step consumes
numbers from it, so a run that skips the HMM would shift every later split. The MLP-only mode replays those draws
(`replay_hmm_rng`), and a dataset arm that is switched off does the same. The per-iteration sizes the runs print
(`grep -E "full train=|balanced train=" <log> | md5sum`) have the same fingerprint, `c4f43c0c9a04cd2e6d46ada7dc7e7327`,
in the original run, in both re-runs above and in the validation run below, for all 32 iterations (16 cluster-disjoint and
16 random), so the shared test sets are identical.

## How large the run-to-run noise is

The torch random state is consumed across the MLP variants, so a re-run is a replicate rather than a bit-identical copy.
For MLP_ProstT5 with the same configuration, the cluster-split mean MCC was 0.928 (original run), 0.926 (job 300283, the
committed rows) and 0.933 (job 300285, below); per-iteration differences between runs average about 0.008 and reach
0.04. The ranking of the top five MLPs (0.919–0.928 here) therefore changes between runs, while their statistical tie does
not depend on it.

Validation run: a single MLP-only run of all six variants with the final script defaults
(`run_bakeoff_mlp_only.sh`, SLURM job 300285; 11 minutes on one GPU) gave cluster-split mean MCC within 0.007 of the
committed values:

| MLP | committed | job 300285 |
|---|---|---|
| ProstT5 | 0.926 | 0.933 |
| SaProt | 0.927 | 0.931 |
| ProtT5+ESM-IF | 0.928 | 0.925 |
| ProtT5 | 0.919 | 0.923 |
| ProstT5+ESM-IF | 0.920 | 0.921 |
| ESM-IF | 0.859 | 0.863 |

(random split: all within 0.004).

## Reproducing

```bash
cd 2_training
sbatch run_bakeoff_mlp_only.sh      # MLP rows only, same splits; then swap them in with merge_bakeoff_rows.py
sbatch run_comparison.sh            # everything from scratch (hours; needs mmseqs, foldseek, mafft, hmmer on PATH)
```
