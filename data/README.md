# Data

Small reference and result data files included in the repository (all <10 MB).

## Files

| File | Size | Description |
|------|------|-------------|
| `classifier_results_all_datasets.csv` | 88 KB | Per-iteration classifier performance metrics (AUPRC, AUROC, MCC, F1, Accuracy) across all datasets and clustering thresholds. |
| `classifier_results_balanced_datasets.csv` | 13 KB | Classifier performance on class-balanced evaluation datasets. |
| `model_comparison_summary.csv` | <1 KB | Summary table comparing MLP vs. XGBoost performance. |
| `reference_proteomes.csv` | 222 KB | Metadata for 386 reference proteomes (organism, taxonomy, proteome ID). |

## Large Files (Not Included)

The following files are excluded via `.gitignore` due to size:

- `SAM_swiss_substrate_types.csv` (14 MB) — UniProtKB/Swiss-Prot SAM MTase annotations with substrate labels
- `*.h5` — ProtT5 embedding files (often >1 GB)
- `*.fasta` — Protein sequence files
- `predictions_G_df*.csv` — Genome-scale prediction output files

These files can be regenerated using the scripts in `embeddings/` and `predictions/`.
