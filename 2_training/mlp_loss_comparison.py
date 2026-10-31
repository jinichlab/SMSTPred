#!/usr/bin/env python3
"""
MLP loss-function & training-epoch comparison (consolidated, parametrized).

This single script performs the loss/epoch tuning stage for any of the three
training sets. It replaces three byte-for-byte-identical scripts that differed
only in a small config block:
  - SAM_mlp_weighted_bce_comparison.py            (DATASET="SAM")
  - SAM_balanced_mlp_weighted_bce_comparison.py   (DATASET="SAM_balanced")
  - UP_mlp_weighted_bce_comparison.py             (DATASET="UP")

What it does
------------
Under cluster-disjoint cross-validation (3 iterations), the MLP
(1024 -> 512 -> 256 -> 4) is trained up to 20,000 epochs and evaluated at
1,000-epoch checkpoints, for each of four candidate per-label-uniform losses:
  'bce', 'weighted_bce', 'sqrt_weighted_bce', 'focal'.

The payoff is the PER-LABEL CSV: it shows which loss wins for each label
(Metabolite / Protein / RNA / DNA), which is what motivates the per-label loss
combination used by the deployment pipeline (e.g. SAMImbalancedLoss). The
`best_*_model.pth` saved here is just the best uniform-loss checkpoint, NOT the
deployed model.

Usage
-----
  python mlp_loss_comparison.py              # run all three datasets
  python mlp_loss_comparison.py SAM          # run a single dataset
  python mlp_loss_comparison.py UP

Outputs (per dataset, written to ./<DATASET>_loss_comparison_output/)
  - <DATASET>_mlp_loss_comparison_summary.csv      mean +/- SD per loss/epoch
  - <DATASET>_mlp_loss_comparison_detailed.csv     per-iteration metrics
  - <DATASET>_mlp_loss_comparison_per_label.csv    per-label breakdown
  - <DATASET>_mlp_loss_comparison_training_curves.png
  - <DATASET>_mlp_loss_comparison_metrics.png
  - best_<DATASET>_mlp_loss_comparison_model.pth
"""

import os
import sys
import pickle
import random
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import (
    matthews_corrcoef,
    precision_recall_curve,
    auc,
    roc_auc_score,
    f1_score,
    accuracy_score,
)
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import matplotlib.pyplot as plt
import psutil

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)


# ============================================================================
# Configuration  --  data directory (see data/README.md; override with $SMSTPRED_DATA)
# ============================================================================
DATA_ROOT = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[1] / "data"))
FILES_DIR = DATA_ROOT

# Per-dataset config. The three training sets differ only in these values.
DATASETS = {
    "SAM": {
        "csv": FILES_DIR / "labels_SAM_full.csv",
        "h5": FILES_DIR / "embeddings_SAM_ProtT5.h5",
        "threshold": 0.6,
    },
    "SAM_balanced": {
        "csv": FILES_DIR / "labels_SAM_balanced.csv",
        "h5": FILES_DIR / "embeddings_SAM_ProtT5.h5",
        "threshold": 0.6,
    },
    "UP": {
        "csv": FILES_DIR / "labels_UP.csv",
        "h5": FILES_DIR / "embeddings_UP_ProtT5.h5",
        "threshold": 0.5,
    },
}

LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
LOSS_TYPES = ["bce", "weighted_bce", "sqrt_weighted_bce", "focal"]
MAX_EPOCHS = 20000
EPOCH_CHECKPOINTS = list(range(1000, 21000, 1000))
LEARNING_RATE = 0.001
N_ITERATIONS = 3


# ============================================================================
# Data Loading
# ============================================================================

def get_embeds(prott5_file, df_file):
    """Load embeddings from H5 file and merge with dataframe."""
    data_map = {}
    with h5py.File(prott5_file, "r") as file:
        for sequence_id, embedding in file.items():
            data_map[sequence_id] = np.array(embedding)

    prott5_df = df_file.copy()
    prott5_df["Embeddings"] = prott5_df["Entry"].map(data_map)
    prott5_df["Embeddings"] = prott5_df["Embeddings"].apply(
        lambda x: x.tolist() if x is not None else None
    )
    prott5_df = prott5_df.dropna(subset=["Embeddings"])
    return prott5_df


# ============================================================================
# Clustering Functions
# ============================================================================

def estimate_memory_usage(n_samples):
    return (n_samples ** 2 * 4) / (1024 ** 3)


def batched_cosine_similarity(embeddings, batch_size=5000, save_path="cosine_sim_matrix.dat"):
    """Compute cosine similarity in batches & store in disk-backed array."""
    n = embeddings.shape[0]
    if os.path.exists(save_path):
        print(f"Loading existing similarity matrix from {save_path}")
        return np.memmap(save_path, dtype=np.float32, mode="r", shape=(n, n))
    sim_matrix = np.memmap(save_path, dtype=np.float32, mode="w+", shape=(n, n))
    for i in tqdm(range(0, n, batch_size), desc="Computing Cosine Similarity", leave=True):
        end_i = min(i + batch_size, n)
        for j in range(0, n, batch_size):
            end_j = min(j + batch_size, n)
            sim_matrix[i:end_i, j:end_j] = cosine_similarity(
                embeddings[i:end_i], embeddings[j:end_j]
            )
    return sim_matrix


def on_demand_cosine_similarity(embeddings, threshold):
    """Greedy clustering with on-demand similarity calculation."""
    n = len(embeddings)
    embedding_norms = np.linalg.norm(embeddings, axis=1)
    sorted_indices = np.argsort(-embedding_norms)

    clusters = {}
    cluster_labels = -1 * np.ones(n, dtype=int)
    cluster_id = 0
    for i in tqdm(sorted_indices, desc="On-Demand Greedy Clustering", leave=True):
        if cluster_labels[i] != -1:
            continue
        assigned = False
        for rep in clusters.keys():
            rep_index = clusters[rep][0]
            if cosine_similarity([embeddings[i]], [embeddings[rep_index]])[0, 0] >= threshold:
                cluster_labels[i] = rep
                clusters[rep].append(i)
                assigned = True
                break
        if not assigned:
            clusters[cluster_id] = [i]
            cluster_labels[i] = cluster_id
            cluster_id += 1
    return clusters, cluster_labels


def greedy_clustering_from_matrix(sim_matrix, threshold):
    """Greedy clustering using a precomputed cosine similarity matrix."""
    n = len(sim_matrix)
    embedding_norms = np.linalg.norm(sim_matrix, axis=1)
    sorted_indices = np.argsort(-embedding_norms)

    clusters = {}
    cluster_labels = -1 * np.ones(n, dtype=int)
    cluster_id = 0
    for i in tqdm(sorted_indices, desc="Greedy Clustering (Precomputed)", leave=True):
        if cluster_labels[i] != -1:
            continue
        assigned = False
        for rep in clusters.keys():
            if sim_matrix[i, clusters[rep][0]] >= threshold:
                cluster_labels[i] = rep
                clusters[rep].append(i)
                assigned = True
                break
        if not assigned:
            clusters[cluster_id] = [i]
            cluster_labels[i] = cluster_id
            cluster_id += 1
    return clusters, cluster_labels


def dynamic_clustering(embeddings, threshold, sim_matrix_path="cosine_sim_matrix.dat", cluster_path=None):
    """Store full similarity matrix if memory allows, else compute on demand.
    Reuses a cached cluster pickle if present (e.g. from the clustering stage)."""
    if cluster_path and os.path.exists(cluster_path):
        print(f"\nLoading cached clusters from {cluster_path}")
        with open(cluster_path, "rb") as f:
            cached = pickle.load(f)
        print(f"  Loaded {len(cached['clusters'])} clusters")
        return cached["clusters"], cached["cluster_labels"]

    n_samples = len(embeddings)
    available_memory_gb = psutil.virtual_memory().available / (1024 ** 3)
    required_memory_gb = estimate_memory_usage(n_samples)
    print(f"\nDataset Size: {n_samples} samples")
    print(f"Memory Needed for Full Matrix: {required_memory_gb:.2f} GB")
    print(f"Available System Memory: {available_memory_gb:.2f} GB")

    if required_memory_gb < available_memory_gb * 0.5:
        print("Using Precomputed Cosine Similarity Matrix")
        sim_matrix = batched_cosine_similarity(embeddings, 5000, save_path=sim_matrix_path)
        clusters, cluster_labels = greedy_clustering_from_matrix(sim_matrix, threshold)
    else:
        print("Memory Constraint! Switching to On-Demand Similarity Calculation")
        clusters, cluster_labels = on_demand_cosine_similarity(embeddings, threshold)

    if cluster_path:
        with open(cluster_path, "wb") as f:
            pickle.dump({"clusters": clusters, "cluster_labels": cluster_labels}, f)
        print(f"Clusters saved to {cluster_path}")
    return clusters, cluster_labels


# ============================================================================
# Splitting Functions
# ============================================================================

def split_clusters_train_val(clusters, y_labels, train_ratio=0.9, max_attempts=1000):
    """Split clusters into train (90%) / val (10%), ensuring every label has at
    least one positive in validation."""
    cluster_ids = list(clusters.keys())
    n_labels = y_labels.shape[1]
    for attempt in range(max_attempts):
        random.shuffle(cluster_ids)
        train_size = int(train_ratio * len(cluster_ids))
        train_clusters = {k: clusters[k] for k in cluster_ids[:train_size]}
        val_clusters = {k: clusters[k] for k in cluster_ids[train_size:]}
        train_indices = [idx for cluster in train_clusters.values() for idx in cluster]
        val_indices = [idx for cluster in val_clusters.values() for idx in cluster]
        y_val = y_labels[val_indices]
        if all(y_val[:, j].sum() > 0 for j in range(n_labels)):
            if attempt > 0:
                print(f"    Found valid split after {attempt + 1} attempts")
            return train_indices, val_indices, train_clusters
    raise ValueError(
        f"Could not find a valid cluster split after {max_attempts} attempts."
    )


def split_train_pool_for_epoch(n_samples, train_ratio=0.9):
    """Randomly split training pool into train (90%) / test (10%) for one epoch."""
    all_indices = list(range(n_samples))
    random.shuffle(all_indices)
    train_size = int(train_ratio * n_samples)
    return all_indices[:train_size], all_indices[train_size:]


# ============================================================================
# Neural Network Model
# ============================================================================

class MultilabelNN(nn.Module):
    def __init__(self, input_size, output_size, use_sigmoid=True):
        super().__init__()
        self.fc1 = nn.Linear(input_size, 512)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, output_size)
        self.use_sigmoid = use_sigmoid

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.relu(self.fc2(x))
        x = self.fc3(x)
        if self.use_sigmoid:
            x = torch.sigmoid(x)
        return x


# ============================================================================
# Metric Computation
# ============================================================================

def compute_auprc(y_true, y_pred_probs):
    auprcs = []
    for idx in range(y_true.shape[1]):
        if len(np.unique(y_true[:, idx])) > 1:
            precision, recall, _ = precision_recall_curve(y_true[:, idx], y_pred_probs[:, idx])
            auprcs.append(auc(recall, precision))
    return np.mean(auprcs) if auprcs else 0.0


def compute_auroc(y_true, y_pred_probs):
    aurocs = []
    for idx in range(y_true.shape[1]):
        if len(np.unique(y_true[:, idx])) > 1:
            aurocs.append(roc_auc_score(y_true[:, idx], y_pred_probs[:, idx]))
    return np.mean(aurocs) if aurocs else 0.0


def compute_mcc(y_true, y_pred):
    mccs = []
    for idx in range(y_true.shape[1]):
        if len(np.unique(y_true[:, idx])) > 1:
            mccs.append(matthews_corrcoef(y_true[:, idx], y_pred[:, idx]))
    return np.mean(mccs) if mccs else 0.0


def compute_f1(y_true, y_pred):
    f1s = []
    for idx in range(y_true.shape[1]):
        if len(np.unique(y_true[:, idx])) > 1:
            f1s.append(f1_score(y_true[:, idx], y_pred[:, idx]))
    return np.mean(f1s) if f1s else 0.0


def compute_accuracy(y_true, y_pred):
    accs = [accuracy_score(y_true[:, idx], y_pred[:, idx]) for idx in range(y_true.shape[1])]
    return np.mean(accs) if accs else 0.0


def compute_all_metrics(y_true, y_pred_probs, threshold=0.5):
    y_pred = (y_pred_probs >= threshold).astype(int)
    return {
        "auprc": compute_auprc(y_true, y_pred_probs),
        "mcc": compute_mcc(y_true, y_pred),
        "auroc": compute_auroc(y_true, y_pred_probs),
        "f1": compute_f1(y_true, y_pred),
        "accuracy": compute_accuracy(y_true, y_pred),
    }


def compute_per_label_metrics(y_true, y_pred_probs, labels, threshold=0.5):
    y_pred = (y_pred_probs >= threshold).astype(int)
    per_label = {}
    for idx, label in enumerate(labels):
        if len(np.unique(y_true[:, idx])) > 1:
            precision, recall, _ = precision_recall_curve(y_true[:, idx], y_pred_probs[:, idx])
            per_label[label] = {
                "auprc": auc(recall, precision),
                "auroc": roc_auc_score(y_true[:, idx], y_pred_probs[:, idx]),
                "mcc": matthews_corrcoef(y_true[:, idx], y_pred[:, idx]),
                "f1": f1_score(y_true[:, idx], y_pred[:, idx]),
                "accuracy": accuracy_score(y_true[:, idx], y_pred[:, idx]),
            }
        else:
            per_label[label] = {k: np.nan for k in ["auprc", "auroc", "mcc", "f1", "accuracy"]}
    return per_label


# ============================================================================
# Imbalance-Aware Loss
# ============================================================================

def compute_pos_weight(y_train):
    """Per-label positive weights (neg/pos ratio) for weighted BCE."""
    label_counts = y_train.sum(axis=0).clip(min=1)
    neg_counts = len(y_train) - label_counts
    return torch.tensor(neg_counts / label_counts, dtype=torch.float32)


class PerLabelUniformLoss(nn.Module):
    """Apply the same loss to each label and sum (no averaging across labels).

    loss_type: 'bce' | 'weighted_bce' | 'sqrt_weighted_bce' | 'focal'
    """
    def __init__(self, loss_type, pos_weight=None, gamma=2.0):
        super().__init__()
        self.loss_type = loss_type
        self.gamma = gamma
        if pos_weight is not None:
            self.register_buffer("pos_weight", pos_weight)
            self.register_buffer("sqrt_pos_weight", torch.sqrt(pos_weight))
        else:
            self.pos_weight = None
            self.sqrt_pos_weight = None

    def forward(self, logits, targets):
        total_loss = 0.0
        for i in range(logits.shape[1]):
            logit_i = logits[:, i]
            target_i = targets[:, i]
            if self.loss_type == "bce":
                total_loss += nn.functional.binary_cross_entropy_with_logits(logit_i, target_i)
            elif self.loss_type == "weighted_bce":
                total_loss += nn.functional.binary_cross_entropy_with_logits(
                    logit_i, target_i, pos_weight=self.pos_weight[i:i + 1])
            elif self.loss_type == "sqrt_weighted_bce":
                total_loss += nn.functional.binary_cross_entropy_with_logits(
                    logit_i, target_i, pos_weight=self.sqrt_pos_weight[i:i + 1])
            elif self.loss_type == "focal":
                bce = nn.functional.binary_cross_entropy_with_logits(
                    logit_i, target_i, reduction="none")
                prob = torch.sigmoid(logit_i)
                p_t = target_i * prob + (1 - target_i) * (1 - prob)
                focal_weight = (1 - p_t) ** self.gamma
                pw = target_i * self.pos_weight[i] + (1 - target_i)
                total_loss += (focal_weight * bce * pw).mean()
        return total_loss


# ============================================================================
# MLP Training with Checkpoint Evaluation
# ============================================================================

def train_mlp_with_checkpoints(X_train_full, y_train_full, X_val, y_val,
                               labels, iteration,
                               max_epochs=20000, epoch_checkpoints=None, lr=0.001,
                               loss_type="weighted_bce"):
    if epoch_checkpoints is None:
        epoch_checkpoints = list(range(1000, 21000, 1000))

    print(f"\n  Training MLP (Iteration {iteration + 1}, loss={loss_type}) up to {max_epochs} epochs...")
    torch.set_num_threads(20)

    model = MultilabelNN(input_size=X_train_full.shape[1], output_size=len(labels),
                         use_sigmoid=False)
    pos_weight = compute_pos_weight(y_train_full)
    criterion = PerLabelUniformLoss(loss_type, pos_weight=pos_weight, gamma=2.0)
    print(f"  pos_weight (raw):  {pos_weight.numpy().round(2)}")
    optimizer = optim.Adam(model.parameters(), lr=lr)

    train_losses, test_losses, test_auprcs = [], [], []
    X_val_tensor = torch.tensor(X_val, dtype=torch.float32)
    checkpoint_results = {}
    checkpoint_set = set(epoch_checkpoints)
    n_train = X_train_full.shape[0]

    for epoch in tqdm(range(max_epochs), desc=f"  MLP Training (Iter {iteration + 1})", leave=False):
        train_indices, test_indices = split_train_pool_for_epoch(n_train, train_ratio=0.9)

        X_train_tensor = torch.tensor(X_train_full[train_indices], dtype=torch.float32)
        y_train_tensor = torch.tensor(y_train_full[train_indices], dtype=torch.float32)
        X_test_tensor = torch.tensor(X_train_full[test_indices], dtype=torch.float32)
        y_test = y_train_full[test_indices]
        y_test_tensor = torch.tensor(y_test, dtype=torch.float32)

        model.train()
        optimizer.zero_grad()
        outputs_train = model(X_train_tensor)
        loss_train = criterion(outputs_train, y_train_tensor)
        loss_train.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            outputs_test = model(X_test_tensor)
            loss_test = criterion(outputs_test, y_test_tensor).item()
            y_test_probs = torch.sigmoid(outputs_test).numpy()
            test_auprc = compute_auprc(y_test, y_test_probs)

        train_losses.append(loss_train.item())
        test_losses.append(loss_test)
        test_auprcs.append(test_auprc)

        current_epoch = epoch + 1
        if current_epoch in checkpoint_set:
            model.eval()
            with torch.no_grad():
                y_val_probs = torch.sigmoid(model(X_val_tensor)).numpy()
            metrics = compute_all_metrics(y_val, y_val_probs)
            per_label = compute_per_label_metrics(y_val, y_val_probs, labels)
            checkpoint_results[current_epoch] = {
                "metrics": metrics,
                "per_label": per_label,
                "model_state": {k: v.clone() for k, v in model.state_dict().items()},
            }
            print(f"\n    Epoch {current_epoch}: AUPRC={metrics['auprc']:.4f}  "
                  f"AUROC={metrics['auroc']:.4f}  MCC={metrics['mcc']:.4f}  "
                  f"F1={metrics['f1']:.4f}  Acc={metrics['accuracy']:.4f}")

    return model, checkpoint_results, train_losses, test_losses, test_auprcs


# ============================================================================
# Per-dataset run
# ============================================================================

# Reuse the cluster assignments produced by the Stage 1 clustering scripts when
# available (same training set, same threshold -> identical partition), matching
# the original workflow where the tuning step shared the sweep's cluster pickle.
# SAM_balanced has no Stage 1 equivalent (it is a subset), so it always
# reclusters locally.
STAGE1_DIR = Path(__file__).resolve().parent.parent / "1_clustering"
STAGE1_PICKLES = {
    "SAM": STAGE1_DIR / "sam_clustering_output" / "SAM_clusters_t0.6.pkl",
    "UP": STAGE1_DIR / "up_clustering_output" / "UP_clusters_t0.5.pkl",
}


def run_dataset(name, cfg):
    out_dir = Path(__file__).resolve().parent / f"{name}_loss_comparison_output"
    out_dir.mkdir(parents=True, exist_ok=True)
    threshold = cfg["threshold"]

    print("=" * 70)
    print(f"{name} MLP loss comparison: BCE vs Weighted vs Sqrt-Dampened vs Focal")
    print(f"Checkpoints: {EPOCH_CHECKPOINTS}")
    print("=" * 70)

    print("\nLoading training dataset...")
    df = pd.read_csv(cfg["csv"])
    df = get_embeds(str(cfg["h5"]), df)
    print(f"Dataset shape after loading embeddings: {df.shape}")

    embeddings = np.vstack(df["Embeddings"].values).astype(np.float32)
    y_all = df[LABELS].values.astype(np.float32)

    # Prefer a Stage 1 cluster pickle if present (avoids recomputing the heavy
    # UP similarity matrix + greedy clustering a second time); else cluster locally.
    stage1 = STAGE1_PICKLES.get(name)
    if stage1 and stage1.exists():
        cluster_path = str(stage1)
        print(f"\nReusing Stage 1 clusters: {cluster_path}")
    else:
        cluster_path = str(out_dir / f"{name}_clusters_t{threshold}.pkl")
        print(f"\nClustering at threshold {threshold}...")
    clusters, _ = dynamic_clustering(
        embeddings, threshold,
        sim_matrix_path=str(out_dir / f"{name}_cosine_sim_matrix.dat"),
        cluster_path=cluster_path,
    )
    print(f"Number of clusters: {len(clusters)}")

    all_results, per_label_results = [], []
    best_model_state, best_auprc, best_epoch, best_iteration, best_loss_type = None, -1, -1, -1, None
    all_training_curves = {}

    for iteration in range(N_ITERATIONS):
        print("\n" + "=" * 70)
        print(f"  ITERATION {iteration + 1}/{N_ITERATIONS}")
        print("=" * 70)

        train_indices, val_indices, _ = split_clusters_train_val(clusters, y_all, train_ratio=0.9)
        print(f"  Train samples: {len(train_indices)}, Validation samples: {len(val_indices)}")

        df_train = df.iloc[train_indices].reset_index(drop=True)
        X_train = np.vstack(df_train["Embeddings"].values).astype(np.float32)
        y_train = df_train[LABELS].values.astype(np.float32)
        df_val = df.iloc[val_indices].reset_index(drop=True)
        X_val = np.vstack(df_val["Embeddings"].values).astype(np.float32)
        y_val = df_val[LABELS].values.astype(np.float32)

        for loss_type in LOSS_TYPES:
            model, checkpoint_results, train_losses, test_losses, test_auprcs = \
                train_mlp_with_checkpoints(
                    X_train, y_train, X_val, y_val, LABELS, iteration,
                    max_epochs=MAX_EPOCHS, epoch_checkpoints=EPOCH_CHECKPOINTS,
                    lr=LEARNING_RATE, loss_type=loss_type,
                )
            all_training_curves[(iteration, loss_type)] = (train_losses, test_losses, test_auprcs)

            for epoch, ckpt in checkpoint_results.items():
                metrics = ckpt["metrics"]
                all_results.append({"epochs": epoch, "iteration": iteration + 1,
                                    "loss_type": loss_type, **metrics})
                for label_name, label_metrics in ckpt["per_label"].items():
                    per_label_results.append({"epochs": epoch, "iteration": iteration + 1,
                                              "loss_type": loss_type, "label": label_name,
                                              **label_metrics})
                if metrics["auprc"] > best_auprc:
                    best_auprc, best_model_state = metrics["auprc"], ckpt["model_state"]
                    best_epoch, best_iteration, best_loss_type = epoch, iteration + 1, loss_type

    # ---- Summaries ----
    results_df = pd.DataFrame(all_results)
    metrics_names = ["auprc", "auroc", "mcc", "f1", "accuracy"]
    summary_data = []
    for loss_type in LOSS_TYPES:
        subset_loss = results_df[results_df["loss_type"] == loss_type]
        for epoch in EPOCH_CHECKPOINTS:
            subset = subset_loss[subset_loss["epochs"] == epoch]
            for metric in metrics_names:
                vals = subset[metric].values
                summary_data.append({"LossType": loss_type, "Epochs": epoch,
                                     "Metric": metric.upper(),
                                     "Mean": np.mean(vals), "Std": np.std(vals)})

    pd.DataFrame(summary_data).to_csv(out_dir / f"{name}_mlp_loss_comparison_summary.csv", index=False)
    results_df.to_csv(out_dir / f"{name}_mlp_loss_comparison_detailed.csv", index=False)
    pd.DataFrame(per_label_results).to_csv(out_dir / f"{name}_mlp_loss_comparison_per_label.csv", index=False)
    print(f"\nSummary / detailed / per-label CSVs written to {out_dir}")

    # ---- Training-curve figure ----
    loss_colors = {"bce": "forestgreen", "weighted_bce": "darkorange",
                   "sqrt_weighted_bce": "mediumorchid", "focal": "steelblue"}
    n_rows = len(LOSS_TYPES) * 2
    fig, axes = plt.subplots(n_rows, N_ITERATIONS, figsize=(6 * N_ITERATIONS, 5 * n_rows))
    for row_block, loss_type in enumerate(LOSS_TYPES):
        for i in range(N_ITERATIONS):
            train_losses, test_losses, test_auprcs = all_training_curves[(i, loss_type)]
            color = loss_colors[loss_type]
            ax1 = axes[row_block * 2, i]
            ax1.plot(range(1, MAX_EPOCHS + 1), train_losses, label="Train Loss", color=color, alpha=0.7)
            ax1.plot(range(1, MAX_EPOCHS + 1), test_losses, label="Test Loss", color=color, alpha=0.4, linestyle="--")
            for ep in EPOCH_CHECKPOINTS:
                ax1.axvline(x=ep, color="gray", linestyle="--", alpha=0.5)
            ax1.set_xlabel("Epoch"); ax1.set_ylabel("Loss")
            ax1.set_title(f"{loss_type.upper()} - Iter {i + 1}: Loss"); ax1.legend(); ax1.grid(True, alpha=0.3)
            ax2 = axes[row_block * 2 + 1, i]
            ax2.plot(range(1, MAX_EPOCHS + 1), test_auprcs, label="Test AUPRC", color=color, alpha=0.7)
            for ep in EPOCH_CHECKPOINTS:
                ax2.axvline(x=ep, color="gray", linestyle="--", alpha=0.5)
                sub = results_df[(results_df["epochs"] == ep) & (results_df["iteration"] == i + 1) &
                                 (results_df["loss_type"] == loss_type)]
                if len(sub) > 0:
                    ax2.plot(ep, sub["auprc"].values[0], "ro", markersize=8)
            ax2.set_xlabel("Epoch"); ax2.set_ylabel("AUPRC")
            ax2.set_title(f"{loss_type.upper()} - Iter {i + 1}: AUPRC"); ax2.legend(); ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / f"{name}_mlp_loss_comparison_training_curves.png", dpi=150)
    plt.close()

    # ---- Metric-comparison figure ----
    fig, axes = plt.subplots(1, len(metrics_names), figsize=(4 * len(metrics_names), 5))
    for idx, metric in enumerate(metrics_names):
        ax = axes[idx]
        for loss_type in LOSS_TYPES:
            subset_loss = results_df[results_df["loss_type"] == loss_type]
            means = [np.mean(subset_loss[subset_loss["epochs"] == ep][metric].values) for ep in EPOCH_CHECKPOINTS]
            stds = [np.std(subset_loss[subset_loss["epochs"] == ep][metric].values) for ep in EPOCH_CHECKPOINTS]
            ax.errorbar(EPOCH_CHECKPOINTS, means, yerr=stds, marker="o", capsize=5,
                        label=loss_type.upper(), color=loss_colors[loss_type])
        ax.set_xlabel("Epochs"); ax.set_ylabel(metric.upper()); ax.set_title(metric.upper())
        ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / f"{name}_mlp_loss_comparison_metrics.png", dpi=150)
    plt.close()

    torch.save(best_model_state, out_dir / f"best_{name}_mlp_loss_comparison_model.pth")
    print(f"\nBest model: iteration {best_iteration}, epoch {best_epoch}, loss={best_loss_type} "
          f"(AUPRC: {best_auprc:.4f}) -> best_{name}_mlp_loss_comparison_model.pth")
    print("=" * 70 + f"\n{name} COMPLETE\n" + "=" * 70)


# ============================================================================
# Main
# ============================================================================

def main():
    if len(sys.argv) > 1:
        requested = sys.argv[1]
        if requested not in DATASETS:
            sys.exit(f"Unknown dataset '{requested}'. Choose from: {', '.join(DATASETS)}")
        to_run = [requested]
    else:
        to_run = list(DATASETS.keys())

    for name in to_run:
        run_dataset(name, DATASETS[name])


if __name__ == "__main__":
    main()
