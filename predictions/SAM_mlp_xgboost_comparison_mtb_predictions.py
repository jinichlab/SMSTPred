import os
import gzip
import random
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
from xgboost import XGBClassifier
from sklearn.multioutput import MultiOutputClassifier
from tqdm import tqdm
import matplotlib.pyplot as plt
import psutil

# Set random seeds for reproducibility
random.seed(42)
np.random.seed(42)
torch.manual_seed(42)


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
    prott5_df['Embeddings'] = prott5_df['Entry'].map(data_map)
    prott5_df['Embeddings'] = prott5_df['Embeddings'].apply(
        lambda x: x.tolist() if x is not None else None
    )
    prott5_df = prott5_df.dropna(subset=['Embeddings'])
    return prott5_df


# ============================================================================
# Clustering Functions
# ============================================================================

def estimate_memory_usage(n_samples):
    """Estimate memory required for a full cosine similarity matrix."""
    memory_needed_gb = (n_samples ** 2 * 4) / (1024 ** 3)
    return memory_needed_gb


def batched_cosine_similarity(embeddings, batch_size=5000, save_path="cosine_sim_matrix.dat"):
    """Compute cosine similarity in batches & store in disk-backed array."""
    n = embeddings.shape[0]
    sim_matrix = np.memmap(save_path, dtype=np.float32, mode='w+', shape=(n, n))

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

        clusters[cluster_id] = [i]
        cluster_labels[i] = cluster_id

        for rep in clusters.keys():
            rep_index = clusters[rep][0]
            sim = cosine_similarity([embeddings[i]], [embeddings[rep_index]])[0, 0]

            if sim >= threshold:
                cluster_labels[i] = rep
                clusters[rep].append(i)
                break

        if cluster_labels[i] == cluster_id:
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

        clusters[cluster_id] = [i]
        cluster_labels[i] = cluster_id

        for rep in clusters.keys():
            rep_index = clusters[rep][0]
            sim = sim_matrix[i, rep_index]

            if sim >= threshold:
                cluster_labels[i] = rep
                clusters[rep].append(i)
                break

        if cluster_labels[i] == cluster_id:
            cluster_id += 1

    return clusters, cluster_labels


def dynamic_clustering(embeddings, threshold):
    """Automatically determines whether to store full similarity matrix or compute on demand."""
    n_samples = len(embeddings)
    available_memory_gb = psutil.virtual_memory().available / (1024 ** 3)
    required_memory_gb = estimate_memory_usage(n_samples)

    print(f"\nDataset Size: {n_samples} samples")
    print(f"Memory Needed for Full Matrix: {required_memory_gb:.2f} GB")
    print(f"Available System Memory: {available_memory_gb:.2f} GB")

    if required_memory_gb < available_memory_gb * 0.5:
        print("Using Precomputed Cosine Similarity Matrix")
        sim_matrix = batched_cosine_similarity(embeddings, 5000)
        return greedy_clustering_from_matrix(sim_matrix, threshold)
    else:
        print("Memory Constraint! Switching to On-Demand Similarity Calculation")
        return on_demand_cosine_similarity(embeddings, threshold)


# ============================================================================
# Cluster Splitting Functions
# ============================================================================

def split_clusters_train_val(clusters, y_labels, train_ratio=0.9, max_attempts=1000):
    """
    Split clusters into train (90%) and validation (10%).
    Ensures every label has at least one positive representative in the validation set.
    Returns indices for each set.

    Args:
        clusters: Dictionary mapping cluster_id -> list of sample indices
        y_labels: numpy array of shape (n_samples, n_labels) with binary labels
        train_ratio: Fraction of clusters to use for training (default 0.9)
        max_attempts: Maximum number of random shuffles to try (default 1000)

    Raises:
        ValueError: If no valid split found after max_attempts
    """
    cluster_ids = list(clusters.keys())
    n_labels = y_labels.shape[1]

    for attempt in range(max_attempts):
        random.shuffle(cluster_ids)

        train_size = int(train_ratio * len(cluster_ids))
        train_clusters = {k: clusters[k] for k in cluster_ids[:train_size]}
        val_clusters = {k: clusters[k] for k in cluster_ids[train_size:]}

        train_indices = [idx for cluster in train_clusters.values() for idx in cluster]
        val_indices = [idx for cluster in val_clusters.values() for idx in cluster]

        # Check if all labels have at least one positive in validation set
        y_val = y_labels[val_indices]
        all_labels_represented = True
        for label_idx in range(n_labels):
            if y_val[:, label_idx].sum() == 0:
                all_labels_represented = False
                break

        if all_labels_represented:
            if attempt > 0:
                print(f"    Found valid split after {attempt + 1} attempts")
            return train_indices, val_indices, train_clusters

    raise ValueError(f"Could not find a valid cluster split after {max_attempts} attempts. "
                     "Consider using a different train_ratio or checking your data distribution.")


def split_train_pool_for_epoch(train_pool_clusters, train_ratio=0.9):
    """
    Randomly split train_pool clusters into train (90%) and test (10%) for each epoch.
    """
    cluster_ids = list(train_pool_clusters.keys())
    random.shuffle(cluster_ids)

    train_size = int(train_ratio * len(cluster_ids))
    train_clusters = {k: train_pool_clusters[k] for k in cluster_ids[:train_size]}
    test_clusters = {k: train_pool_clusters[k] for k in cluster_ids[train_size:]}

    train_indices = [idx for cluster in train_clusters.values() for idx in cluster]
    test_indices = [idx for cluster in test_clusters.values() for idx in cluster]

    return train_indices, test_indices


# ============================================================================
# Neural Network Model
# ============================================================================

class MultilabelNN(nn.Module):
    def __init__(self, input_size, output_size):
        super(MultilabelNN, self).__init__()
        self.fc1 = nn.Linear(input_size, 512)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, output_size)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.relu(x)
        x = self.fc3(x)
        return self.sigmoid(x)


# ============================================================================
# Metric Computation
# ============================================================================

def compute_auprc(y_true, y_pred_probs):
    """Compute macro-averaged AUPRC across all labels."""
    auprcs = []
    for idx in range(y_true.shape[1]):
        if len(np.unique(y_true[:, idx])) > 1:
            precision, recall, _ = precision_recall_curve(y_true[:, idx], y_pred_probs[:, idx])
            auprcs.append(auc(recall, precision))
    return np.mean(auprcs) if auprcs else 0.0


def compute_auroc(y_true, y_pred_probs):
    """Compute macro-averaged AUROC across all labels."""
    aurocs = []
    for idx in range(y_true.shape[1]):
        if len(np.unique(y_true[:, idx])) > 1:
            aurocs.append(roc_auc_score(y_true[:, idx], y_pred_probs[:, idx]))
    return np.mean(aurocs) if aurocs else 0.0


def compute_mcc(y_true, y_pred):
    """Compute macro-averaged MCC across all labels."""
    mccs = []
    for idx in range(y_true.shape[1]):
        if len(np.unique(y_true[:, idx])) > 1:
            mccs.append(matthews_corrcoef(y_true[:, idx], y_pred[:, idx]))
    return np.mean(mccs) if mccs else 0.0


def compute_f1(y_true, y_pred):
    """Compute macro-averaged F1 across all labels."""
    f1s = []
    for idx in range(y_true.shape[1]):
        if len(np.unique(y_true[:, idx])) > 1:
            f1s.append(f1_score(y_true[:, idx], y_pred[:, idx]))
    return np.mean(f1s) if f1s else 0.0


def compute_accuracy(y_true, y_pred):
    """Compute macro-averaged accuracy across all labels."""
    accs = []
    for idx in range(y_true.shape[1]):
        accs.append(accuracy_score(y_true[:, idx], y_pred[:, idx]))
    return np.mean(accs) if accs else 0.0


def compute_all_metrics(y_true, y_pred_probs, threshold=0.5):
    """Compute all metrics: AUPRC, MCC, AUROC, F1, Accuracy."""
    y_pred = (y_pred_probs >= threshold).astype(int)
    return {
        'auprc': compute_auprc(y_true, y_pred_probs),
        'mcc': compute_mcc(y_true, y_pred),
        'auroc': compute_auroc(y_true, y_pred_probs),
        'f1': compute_f1(y_true, y_pred),
        'accuracy': compute_accuracy(y_true, y_pred),
    }


def compute_per_label_metrics(y_true, y_pred_probs, labels, threshold=0.5):
    """Compute per-label metrics."""
    y_pred = (y_pred_probs >= threshold).astype(int)
    per_label = {}
    for idx, label in enumerate(labels):
        per_label[label] = {}
        if len(np.unique(y_true[:, idx])) > 1:
            precision, recall, _ = precision_recall_curve(y_true[:, idx], y_pred_probs[:, idx])
            per_label[label]['auprc'] = auc(recall, precision)
            per_label[label]['auroc'] = roc_auc_score(y_true[:, idx], y_pred_probs[:, idx])
            per_label[label]['mcc'] = matthews_corrcoef(y_true[:, idx], y_pred[:, idx])
            per_label[label]['f1'] = f1_score(y_true[:, idx], y_pred[:, idx])
            per_label[label]['accuracy'] = accuracy_score(y_true[:, idx], y_pred[:, idx])
        else:
            per_label[label] = {k: np.nan for k in ['auprc', 'auroc', 'mcc', 'f1', 'accuracy']}
    return per_label


# ============================================================================
# Training Functions
# ============================================================================

def train_xgboost(X_train, y_train, X_val, y_val, labels, iteration):
    """Train XGBoost multi-label classifier."""
    print(f"\n  Training XGBoost (Iteration {iteration + 1})...")

    xgb_base = XGBClassifier(
        eval_metric='logloss',
        random_state=42 + iteration,
        n_jobs=-1,
        use_label_encoder=False
    )
    model = MultiOutputClassifier(xgb_base)
    model.fit(X_train, y_train)

    # Predict probabilities
    prob_list = [est.predict_proba(X_val) for est in model.estimators_]
    y_val_probs = np.hstack([p[:, 1].reshape(-1, 1) for p in prob_list])

    # Compute metrics
    metrics = compute_all_metrics(y_val, y_val_probs)
    per_label = compute_per_label_metrics(y_val, y_val_probs, labels)

    return model, y_val_probs, metrics, per_label


def train_mlp(X_train_full, y_train_full, X_val, y_val, train_clusters,
              labels, iteration, epochs=500, lr=0.001):
    """Train MLP with epoch-wise train/test splitting."""
    print(f"\n  Training MLP (Iteration {iteration + 1})...")

    device = torch.device("cpu")
    torch.set_num_threads(20)

    # Initialize model
    model = MultilabelNN(input_size=X_train_full.shape[1], output_size=len(labels))
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)

    # Create a mapping from original indices to train_full indices
    train_full_df_indices = list(range(len(X_train_full)))

    # Training tracking
    train_losses = []
    test_losses = []
    test_auprcs = []

    # Prepare validation tensors
    X_val_tensor = torch.tensor(X_val, dtype=torch.float32)
    y_val_tensor = torch.tensor(y_val, dtype=torch.float32)

    for epoch in tqdm(range(epochs), desc=f"  MLP Training (Iter {iteration + 1})", leave=False):
        # Each epoch: randomly split train_pool clusters into train (90%) and test (10%)
        train_indices, test_indices = split_train_pool_for_epoch(train_clusters, train_ratio=0.9)

        # Prepare train data for this epoch
        X_train = X_train_full[train_indices]
        y_train = y_train_full[train_indices]
        X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
        y_train_tensor = torch.tensor(y_train, dtype=torch.float32)

        # Prepare test data for this epoch
        X_test = X_train_full[test_indices]
        y_test = y_train_full[test_indices]
        X_test_tensor = torch.tensor(X_test, dtype=torch.float32)
        y_test_tensor = torch.tensor(y_test, dtype=torch.float32)

        # Training step
        model.train()
        optimizer.zero_grad()
        outputs_train = model(X_train_tensor)
        loss_train = criterion(outputs_train, y_train_tensor)
        loss_train.backward()
        optimizer.step()

        # Evaluation step (on test split for this epoch)
        model.eval()
        with torch.no_grad():
            outputs_test = model(X_test_tensor)
            loss_test = criterion(outputs_test, y_test_tensor).item()
            y_test_probs = outputs_test.numpy()
            test_auprc = compute_auprc(y_test, y_test_probs)

        # Record metrics
        train_losses.append(loss_train.item())
        test_losses.append(loss_test)
        test_auprcs.append(test_auprc)

    # Final evaluation on validation set
    model.eval()
    with torch.no_grad():
        outputs_val = model(X_val_tensor)
        y_val_probs = outputs_val.numpy()

    # Compute metrics
    metrics = compute_all_metrics(y_val, y_val_probs)
    per_label = compute_per_label_metrics(y_val, y_val_probs, labels)

    return model, y_val_probs, metrics, per_label, train_losses, test_losses, test_auprcs


# ============================================================================
# Main Script
# ============================================================================

def main():
    # Configuration
    THRESHOLD = 0.6
    EPOCHS = 500
    LEARNING_RATE = 0.001
    N_ITERATIONS = 3
    LABELS = ['Other', 'Protein', 'RNA', 'DNA']

    print("=" * 70)
    print("XGBoost vs MLP Comparison with Cluster-Based Cross-Validation")
    print("=" * 70)

    # ========================================================================
    # Load Training Data
    # ========================================================================
    print("\nLoading training dataset...")
    SAM_swiss_df = pd.read_csv('../files/SAM_swiss_enzymes_substrate_types.csv')
    print(f"Original dataset shape: {SAM_swiss_df.shape}")

    SAM_swiss_prott5_file = '../files/SAM_swiss.h5'
    SAM_swiss_df = get_embeds(SAM_swiss_prott5_file, SAM_swiss_df)
    print(f"Dataset shape after loading embeddings: {SAM_swiss_df.shape}")

    # Extract embeddings and labels
    embeddings = np.vstack(SAM_swiss_df['Embeddings'].values).astype(np.float32)
    y_all = SAM_swiss_df[LABELS].values.astype(np.float32)

    # ========================================================================
    # Cluster Once (Keep the same clusters for all iterations)
    # ========================================================================
    print(f"\nClustering at threshold {THRESHOLD}...")
    clusters, cluster_labels = dynamic_clustering(embeddings, THRESHOLD)
    print(f"Number of clusters: {len(clusters)}")

    # ========================================================================
    # Prepare Storage for Results
    # ========================================================================
    xgb_results = []
    mlp_results = []
    best_mlp_model = None
    best_mlp_auprc = -1
    best_mlp_iteration = -1

    all_train_losses = []
    all_test_losses = []
    all_test_auprcs = []

    # ========================================================================
    # Run N Iterations
    # ========================================================================
    for iteration in range(N_ITERATIONS):
        print("\n" + "=" * 70)
        print(f"ITERATION {iteration + 1}/{N_ITERATIONS}")
        print("=" * 70)

        # Re-shuffle and split clusters for this iteration
        # Keeps re-picking until all labels have at least one positive in validation
        train_indices, val_indices, train_clusters = split_clusters_train_val(
            clusters, y_all, train_ratio=0.9
        )
        print(f"Train samples: {len(train_indices)}, Validation samples: {len(val_indices)}")

        # Prepare data
        df_train = SAM_swiss_df.iloc[train_indices].reset_index(drop=True)
        X_train = np.vstack(df_train['Embeddings'].values).astype(np.float32)
        y_train = df_train[LABELS].values.astype(np.float32)

        df_val = SAM_swiss_df.iloc[val_indices].reset_index(drop=True)
        X_val = np.vstack(df_val['Embeddings'].values).astype(np.float32)
        y_val = df_val[LABELS].values.astype(np.float32)

        # Remap train_clusters to use indices within X_train
        # Original: train_clusters[cluster_id] = [original_df_indices]
        # We need: remapped_clusters[cluster_id] = [X_train_indices]
        original_to_new = {orig: new for new, orig in enumerate(train_indices)}
        remapped_clusters = {}
        for cluster_id, orig_indices in train_clusters.items():
            remapped_clusters[cluster_id] = [original_to_new[idx] for idx in orig_indices]

        # --------------------------------------------------------------------
        # Train XGBoost
        # --------------------------------------------------------------------
        xgb_model, xgb_probs, xgb_metrics, xgb_per_label = train_xgboost(
            X_train, y_train, X_val, y_val, LABELS, iteration
        )
        xgb_results.append(xgb_metrics)

        print(f"\n  XGBoost Validation Metrics:")
        print(f"    AUPRC: {xgb_metrics['auprc']:.4f}")
        print(f"    AUROC: {xgb_metrics['auroc']:.4f}")
        print(f"    MCC: {xgb_metrics['mcc']:.4f}")
        print(f"    F1: {xgb_metrics['f1']:.4f}")
        print(f"    Accuracy: {xgb_metrics['accuracy']:.4f}")

        # --------------------------------------------------------------------
        # Train MLP
        # --------------------------------------------------------------------
        mlp_model, mlp_probs, mlp_metrics, mlp_per_label, train_losses, test_losses, test_auprcs = train_mlp(
            X_train, y_train, X_val, y_val, remapped_clusters,
            LABELS, iteration, epochs=EPOCHS, lr=LEARNING_RATE
        )
        mlp_results.append(mlp_metrics)
        all_train_losses.append(train_losses)
        all_test_losses.append(test_losses)
        all_test_auprcs.append(test_auprcs)

        print(f"\n  MLP Validation Metrics:")
        print(f"    AUPRC: {mlp_metrics['auprc']:.4f}")
        print(f"    AUROC: {mlp_metrics['auroc']:.4f}")
        print(f"    MCC: {mlp_metrics['mcc']:.4f}")
        print(f"    F1: {mlp_metrics['f1']:.4f}")
        print(f"    Accuracy: {mlp_metrics['accuracy']:.4f}")

        # Track best MLP model
        if mlp_metrics['auprc'] > best_mlp_auprc:
            best_mlp_auprc = mlp_metrics['auprc']
            best_mlp_model = mlp_model
            best_mlp_iteration = iteration + 1

        # Print per-label metrics
        print(f"\n  Per-Label Metrics (Iteration {iteration + 1}):")
        print(f"  {'Label':<10} {'XGB AUPRC':<12} {'MLP AUPRC':<12} {'XGB MCC':<10} {'MLP MCC':<10}")
        print("  " + "-" * 54)
        for label in LABELS:
            xgb_auprc = xgb_per_label[label].get('auprc', np.nan)
            mlp_auprc = mlp_per_label[label].get('auprc', np.nan)
            xgb_mcc = xgb_per_label[label].get('mcc', np.nan)
            mlp_mcc = mlp_per_label[label].get('mcc', np.nan)
            print(f"  {label:<10} {xgb_auprc:<12.4f} {mlp_auprc:<12.4f} {xgb_mcc:<10.4f} {mlp_mcc:<10.4f}")

    # ========================================================================
    # Summary Statistics
    # ========================================================================
    print("\n" + "=" * 70)
    print("SUMMARY: Mean and Std Across All Iterations")
    print("=" * 70)

    metrics_names = ['auprc', 'auroc', 'mcc', 'f1', 'accuracy']

    print(f"\n{'Metric':<12} {'XGBoost Mean':<14} {'XGBoost Std':<14} {'MLP Mean':<14} {'MLP Std':<14}")
    print("-" * 68)

    summary_data = []
    for metric in metrics_names:
        xgb_vals = [r[metric] for r in xgb_results]
        mlp_vals = [r[metric] for r in mlp_results]
        xgb_mean = np.mean(xgb_vals)
        xgb_std = np.std(xgb_vals)
        mlp_mean = np.mean(mlp_vals)
        mlp_std = np.std(mlp_vals)
        print(f"{metric.upper():<12} {xgb_mean:<14.4f} {xgb_std:<14.4f} {mlp_mean:<14.4f} {mlp_std:<14.4f}")
        summary_data.append({
            'Metric': metric.upper(),
            'XGBoost_Mean': xgb_mean,
            'XGBoost_Std': xgb_std,
            'MLP_Mean': mlp_mean,
            'MLP_Std': mlp_std
        })

    # Save summary to CSV
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv('SAM_model_comparison_summary.csv', index=False)
    print("\nSummary saved to 'SAM_model_comparison_summary.csv'")

    # ========================================================================
    # Plot Training Curves for All Iterations
    # ========================================================================
    fig, axes = plt.subplots(2, N_ITERATIONS, figsize=(6 * N_ITERATIONS, 10))

    for i in range(N_ITERATIONS):
        # Loss plot
        ax1 = axes[0, i] if N_ITERATIONS > 1 else axes[0]
        ax1.plot(range(1, EPOCHS + 1), all_train_losses[i], label='Train Loss', alpha=0.7)
        ax1.plot(range(1, EPOCHS + 1), all_test_losses[i], label='Test Loss', alpha=0.7)
        ax1.set_xlabel('Epoch')
        ax1.set_ylabel('BCE Loss')
        ax1.set_title(f'Iteration {i + 1}: Training/Test Loss')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # AUPRC plot
        ax2 = axes[1, i] if N_ITERATIONS > 1 else axes[1]
        ax2.plot(range(1, EPOCHS + 1), all_test_auprcs[i], label='Test AUPRC', color='green', alpha=0.7)
        ax2.axhline(y=mlp_results[i]['auprc'], color='red', linestyle='--',
                   label=f"Val AUPRC: {mlp_results[i]['auprc']:.4f}")
        ax2.set_xlabel('Epoch')
        ax2.set_ylabel('AUPRC')
        ax2.set_title(f'Iteration {i + 1}: Test AUPRC')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('SAM_mlp_training_curves_all_iterations.png', dpi=150)
    plt.close()
    print("Training curves saved to 'SAM_mlp_training_curves_all_iterations.png'")

    # ========================================================================
    # Load Prediction Data and Make Predictions
    # ========================================================================
    print("\n" + "=" * 70)
    print("PREDICTIONS WITH BEST MLP MODEL")
    print("=" * 70)
    print(f"Best MLP model from iteration {best_mlp_iteration} (AUPRC: {best_mlp_auprc:.4f})")

    # Load prediction dataset
    print("\nLoading prediction dataset...")
    mtb_df = pd.read_excel('../files/uniprotkb_taxonomy_id_83332_AND_IPR0290_2025_06_23.xlsx')
    print(f"mtb_df shape before embeddings: {mtb_df.shape}")

    mtb_prott5_file = '../files/SAM_mtb.h5'
    mtb_df = get_embeds(mtb_prott5_file, mtb_df)
    print(f"G_df shape after loading embeddings: {mtb_df.shape}")

    # Prepare prediction data
    X_pred = np.vstack(mtb_df['Embeddings'].values).astype(np.float32)
    X_pred_tensor = torch.tensor(X_pred, dtype=torch.float32)

    # Make predictions with best MLP model
    best_mlp_model.eval()
    with torch.no_grad():
        pred_probs = best_mlp_model(X_pred_tensor).numpy()

    # Add predictions to dataframe
    for idx, label in enumerate(LABELS):
        mtb_df[f'{label}_Prob'] = pred_probs[:, idx]
        mtb_df[f'{label}_Pred'] = (pred_probs[:, idx] >= 0.5).astype(int)

    # Remove Embeddings column for cleaner output
    output_df = mtb_df.drop(columns=['Embeddings'])

    # Save predictions
    output_df.to_csv('mlp_predictions_mtb_df.csv', index=False)
    print(f"\nPredictions saved to 'mlp_predictions_mtb_df.csv'")
    print(f"Prediction columns: {[f'{l}_Prob' for l in LABELS] + [f'{l}_Pred' for l in LABELS]}")

    # Save best model
    torch.save(best_mlp_model.state_dict(), 'best_SAM_mlp_model.pth')
    print(f"Best MLP model saved to 'best_SAM_mlp_model.pth'")

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
