# Predictions

Genome-scale substrate type prediction scripts for SAM-dependent methyltransferases.

## Scripts

| Script | Description |
|--------|-------------|
| `uniprot_mlp_xgboost_comparison_g_prediction.py` | Trains MLP and XGBoost models on Swiss-Prot SAM MTases and predicts substrate types across genome-scale datasets. Uses cluster-based cross-validation with 3 iterations. Outputs per-label probabilities and binary predictions. |
| `SAM_mlp_xgboost_comparison_mtb_predictions.py` | Same pipeline applied specifically to *Mycobacterium tuberculosis* proteome predictions. |

## Usage

```bash
# Run genome-scale predictions
python uniprot_mlp_xgboost_comparison_g_prediction.py

# Run M. tuberculosis predictions
python SAM_mlp_xgboost_comparison_mtb_predictions.py
```

## Output

- CSV files with per-protein substrate type probabilities (RNA, DNA, Protein, Other)
- Binary classification predictions
- Model performance metrics (AUPRC, AUROC, MCC, F1, Accuracy)
- Trained model weights (.pth)
