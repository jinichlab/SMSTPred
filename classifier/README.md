# Classifier

MLP classifier training and evaluation for SAM MTase substrate type prediction.

## Scripts

| Script | Description |
|--------|-------------|
| `mlp_training_with_validation.py` | Trains the MLP classifier (1024→512→256→4) with cluster-based cross-validation. Uses cosine similarity clustering (threshold 0.7) to prevent homology leakage. 500 epochs, Adam optimizer, BCELoss. |
| `classifier_testing.py` | Evaluates multiple classifiers (XGBoost, Random Forest, Logistic Regression, Neural Network) with dynamic clustering and 5 iterations per threshold. |
| `classifier_testing_balanced.py` | Classifier comparison on class-balanced datasets. |
| `cosine_similarity_histogram.py` | Generates cosine similarity distribution plots between protein embeddings. |
| `cluster_similarity_analysis.py` | Analyzes within-cluster and between-cluster similarity patterns. |
| `substrate_type_testing.sh` | Shell script for running classifier evaluation pipelines. |

## Models

The `models/` directory contains trained model weights:

- `mlp_final_model.pth` (~2.6 MB) — Final trained MLP model

## Metrics

The classifier achieves on held-out validation:

| Metric | Value |
|--------|-------|
| AUPRC  | 0.998 |
| AUROC  | 0.999 |
| MCC    | 0.976 |

## Usage

```bash
# Train the MLP classifier
python mlp_training_with_validation.py

# Run classifier comparisons
python classifier_testing.py
```
