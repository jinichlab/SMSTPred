# Stage 3 — Deploying the classifier across the tree of life

Scores new proteins with the Stage-2 deployment model from their ProstT5 embeddings.

## Script

| Script | Purpose |
|---|---|
| `predict.py` | Loads `2_training/deployment_model/` (6-seed ensemble + frozen per-label thresholds + applicability-domain cutoff) and labels new proteins from one ProstT5 embedding `.h5` (or a directory of them). |
| `run_predict.sh` | SLURM driver. |

## Per-protein output columns

| Column | Meaning |
|---|---|
| `score_<label>` | ensemble-mean probability (mean over the 6 members' sigmoids) |
| `call_<label>` | 1 if `score ≥` the frozen per-label threshold (multi-label) |
| `disagreement` | mean over labels of the 6 members' score spread (SD) — high = shaky |
| `ad_maxsim` | max cosine similarity to the training set (ProstT5 space) |
| `in_domain` | 1 if `ad_maxsim ≥` the AD cutoff (≈0.94); 0 = extrapolating |

Thresholds and the AD flag come from cluster-disjoint **out-of-fold** calibration, so scores
are honest in the remote-homology regime; `disagreement` and `in_domain` flag where to
trust a call.

## Outputs

`predictions_output/` — one CSV per input proteome (~1,259 reference proteomes,
~112k proteins scored). These feed the phylogenetic-trend (Stage 4) and UMAP (Stage 5)
analyses, and define the deployment corpus for the Stage-6 interpretability work.
