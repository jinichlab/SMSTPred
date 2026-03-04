# SAM-Dependent Methyltransferase Substrate Type Classifier

Machine learning pipeline for predicting substrate specificity (RNA, DNA, Protein, Other) of S-adenosylmethionine (SAM)-dependent methyltransferases using ProtT5 protein language model embeddings.

## Overview

This project implements a multi-label MLP classifier trained on ProtT5-XL-U50 embeddings of SAM MTases from UniProtKB/Swiss-Prot. The classifier achieves AUPRC 0.998, AUROC 0.999, and MCC 0.976 on held-out validation data using cluster-based cross-validation to prevent homology leakage.

Predictions are applied genome-wide across 386 reference proteomes (312 bacteria, 13 archaea, 61 eukaryotes), revealing universal substrate trade-offs, lineage-specific specializations, and evolutionary patterns in methyltransferase repertoires.

## Repository Structure

```
├── embeddings/          # ProtT5 embedding generation pipeline
├── classifier/          # Model training & evaluation
├── predictions/         # Genome-scale prediction scripts
├── phylogenetics/       # Phylogenetic tree construction & visualization
├── umap/                # UMAP dimensionality reduction & clustering
├── expression_analysis/ # Human tissue expression analysis
├── solubility/          # Protein solubility prediction wrapper
├── data/                # Small reference/result data files
├── figures/             # Publication figures
├── article/             # Journal article
└── environment.yml      # Conda environment specification
```

## Installation

```bash
# Clone the repository
git clone https://github.com/jinichlab/substrate_type_classifier.git
cd substrate_type_classifier

# Create conda environment
conda env create -f environment.yml
conda activate prott5
```

### Requirements

- Python 3.12+
- PyTorch (with CUDA 12 for GPU acceleration)
- Transformers (for ProtT5-XL-U50)
- scikit-learn, XGBoost
- matplotlib, pandas, numpy, h5py

## Pipeline

1. **Generate embeddings**: Convert protein sequences to 1024-dimensional ProtT5 embeddings (`embeddings/`)
2. **Train classifier**: Train MLP with cluster-based cross-validation (`classifier/`)
3. **Predict**: Apply trained model to reference proteomes (`predictions/`)
4. **Visualize**: Generate phylogenetic trees (`phylogenetics/`) and UMAP projections (`umap/`)
5. **Analyze**: Tissue expression patterns (`expression_analysis/`) and solubility (`solubility/`)

## Model Architecture

- **Input**: 1024-dimensional ProtT5-XL-U50 per-protein embeddings
- **Architecture**: MLP (1024 → 512 → 256 → 4) with ReLU activation and sigmoid output
- **Training**: 500 epochs, Adam optimizer (lr=0.001), BCELoss
- **Validation**: Cluster-based cross-validation (cosine similarity threshold 0.7)
- **Output**: Multi-label probabilities for RNA, DNA, Protein, Other substrates

## Citation

If you use this code, please cite:

> Roethler, N., et al. "Phylogenetic Analysis of SAM-Dependent Methyltransferase Substrate Preferences Across 386 Organisms." (2025).

## License

See repository license for details.
