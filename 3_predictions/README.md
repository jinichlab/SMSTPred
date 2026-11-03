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

## Using the refreshed deployment ensemble

`predict.py` reads the deployment model from `2_training/deployment_model/` (the archived ProstT5 ensemble). `DEPLOY_DIR=<release>/deployment python predict.py --input <ProtT5 embeddings> ...` uses the ProtT5
ensemble of `2_training/deploy_ensemble.py` instead (same file formats; the applicability-domain reference is read from the model folder). The input must then be **ProtT5** vectors: the archived tree-of-life
ProtT5 files (`embeddings_tree_of_life_ProtT5.tar`) give the same calls as the training pipeline's vectors (2_training/README.md, "ProtT5 vectors of the proteomes"). The predictions of Stages 4-7 are not recomputed.

Release 2026_03 with the ProtT5 ensemble (`DEPLOY_DIR=<release>/deployment OUT_DIR=<release>/predictions sbatch run_predict.sh <release>/embeddings/proteomes_ProtT5`, 17 seconds): 115,850 proteins, 56.3 % inside the
applicability domain; calls Metabolite 40.9 %, Protein 17.9 %, RNA 31.7 %, DNA 6.3 %; no label called for 5.0 %, more than one for 1.8 %. By length the DNA call rate rises (2.8 % up to 300 aa, 6.2 %, 13.9 %, 35.8 % at 1,001-2,000 aa)
and the inside-the-domain share falls (54.5 %, 61.1 %, 47.4 %, 38.0 %); 86 % of the long proteins called DNA have DNA-methyltransferase or restriction-modification names (against 1.4 % of the long proteins not called DNA),
so the rise is largely real biology, not a length shortcut; the long proteins are the least certain (median member disagreement 0.030 against 0.005). Against the archived ProstT5 calls (111,087 proteins in both) all four calls
agree for 82.3 %; the ProtT5 ensemble calls more Metabolite (+4,899) and fewer DNA (6,082 against 7,808).
