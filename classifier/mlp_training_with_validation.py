import os
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
)
from sklearn.metrics.pairwise import cosine_similarity
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

def split_clusters_three_way(clusters, train_pool_ratio=0.9):
    """
    Split clusters into train_pool (90%) and validation (10%).
    Returns indices for each set.
    """
    cluster_ids = list(clusters.keys())
    random.shuffle(cluster_ids)

    train_pool_size = int(train_pool_ratio * len(cluster_ids))
    train_pool_clusters = {k: clusters[k] for k in cluster_ids[:train_pool_size]}
    val_clusters = {k: clusters[k] for k in cluster_ids[train_pool_size:]}

    train_pool_indices = [idx for cluster in train_pool_clusters.values() for idx in cluster]
    val_indices = [idx for cluster in val_clusters.values() for idx in cluster]

    return train_pool_indices, val_indices, train_pool_clusters


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


def compute_mcc(y_true, y_pred):
    """Compute macro-averaged MCC across all labels."""
    mccs = []
    for idx in range(y_true.shape[1]):
        if len(np.unique(y_true[:, idx])) > 1:
            mccs.append(matthews_corrcoef(y_true[:, idx], y_pred[:, idx]))
    return np.mean(mccs) if mccs else 0.0


# ============================================================================
# Main Training Script
# ============================================================================

def main():
    # Configuration
    THRESHOLD = 0.7
    EPOCHS = 500
    LEARNING_RATE = 0.001
    LABELS = ['Other', 'Protein', 'RNA', 'DNA']

    # Set device
    torch.set_num_threads(20)
    device = torch.device("cpu")

    # Load data
    print("Loading dataset...")
    UP_swiss_df = pd.read_csv('../files/UP_swiss_enzymes_substrate_types.csv')
    print(f"Original dataset shape: {UP_swiss_df.shape}")

    UP_swiss_prott5_file = '../files/per-protein.h5'
    UP_swiss_df = get_embeds(UP_swiss_prott5_file, UP_swiss_df)
    print(f"Dataset shape after loading embeddings: {UP_swiss_df.shape}")

    # Extract embeddings
    embeddings = np.vstack(UP_swiss_df['Embeddings'].values).astype(np.float32)

    # Cluster at 0.7 threshold
    print(f"\nClustering at threshold {THRESHOLD}...")
    clusters, cluster_labels = dynamic_clustering(embeddings, THRESHOLD)
    print(f"Number of clusters: {len(clusters)}")

    # Split: 90% train_pool, 10% validation (held out)
    print("\nSplitting clusters into train_pool (90%) and validation (10%)...")
    train_pool_indices, val_indices, train_pool_clusters = split_clusters_three_way(
        clusters, train_pool_ratio=0.9
    )
    print(f"Train pool samples: {len(train_pool_indices)}")
    print(f"Validation samples: {len(val_indices)}")

    # Prepare validation data (held out until the end)
    df_val = UP_swiss_df.iloc[val_indices].reset_index(drop=True)
    X_val = np.vstack(df_val['Embeddings'].values).astype(np.float32)
    y_val = df_val[LABELS].values.astype(np.float32)
    X_val_tensor = torch.tensor(X_val, dtype=torch.float32)
    y_val_tensor = torch.tensor(y_val, dtype=torch.float32)

    # Initialize model
    model = MultilabelNN(input_size=1024, output_size=len(LABELS))
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)

    # Training tracking
    train_losses = []
    test_losses = []
    test_auprcs = []

    print(f"\nStarting training for {EPOCHS} epochs...")
    for epoch in tqdm(range(EPOCHS), desc="Training"):
        # Each epoch: randomly split train_pool into train (90%) and test (10%)
        train_indices, test_indices = split_train_pool_for_epoch(
            train_pool_clusters, train_ratio=0.9
        )

        # Prepare train data for this epoch
        df_train = UP_swiss_df.iloc[train_indices].reset_index(drop=True)
        X_train = np.vstack(df_train['Embeddings'].values).astype(np.float32)
        y_train = df_train[LABELS].values.astype(np.float32)
        X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
        y_train_tensor = torch.tensor(y_train, dtype=torch.float32)

        # Prepare test data for this epoch
        df_test = UP_swiss_df.iloc[test_indices].reset_index(drop=True)
        X_test = np.vstack(df_test['Embeddings'].values).astype(np.float32)
        y_test = df_test[LABELS].values.astype(np.float32)
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

        # Print progress every 50 epochs
        if (epoch + 1) % 50 == 0:
            print(f"\nEpoch {epoch + 1}/{EPOCHS} | "
                  f"Train Loss: {loss_train.item():.4f} | "
                  f"Test Loss: {loss_test:.4f} | "
                  f"Test AUPRC: {test_auprc:.4f}")

    # ========================================================================
    # Final Evaluation on Held-Out Validation Set
    # ========================================================================
    print("\n" + "=" * 60)
    print("Final Evaluation on Held-Out Validation Set")
    print("=" * 60)

    model.eval()
    with torch.no_grad():
        outputs_val = model(X_val_tensor)
        y_val_probs = outputs_val.numpy()
        y_val_pred = (y_val_probs >= 0.5).astype(int)

    # Compute metrics
    val_auprc = compute_auprc(y_val, y_val_probs)
    val_mcc = compute_mcc(y_val, y_val_pred)

    print(f"\nValidation AUPRC: {val_auprc:.4f}")
    print(f"Validation MCC: {val_mcc:.4f}")

    # Per-label metrics
    print("\nPer-label metrics:")
    for idx, label in enumerate(LABELS):
        if len(np.unique(y_val[:, idx])) > 1:
            precision, recall, _ = precision_recall_curve(y_val[:, idx], y_val_probs[:, idx])
            label_auprc = auc(recall, precision)
            label_mcc = matthews_corrcoef(y_val[:, idx], y_val_pred[:, idx])
            print(f"  {label}: AUPRC={label_auprc:.4f}, MCC={label_mcc:.4f}")
        else:
            print(f"  {label}: Only one class present in validation set")

    # ========================================================================
    # Plotting
    # ========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Loss across epochs
    axes[0].plot(range(1, EPOCHS + 1), train_losses, label='Train Loss', alpha=0.7)
    axes[0].plot(range(1, EPOCHS + 1), test_losses, label='Test Loss', alpha=0.7)
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('BCE Loss')
    axes[0].set_title('Training and Test Loss Across Epochs')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Plot 2: AUPRC across epochs
    axes[1].plot(range(1, EPOCHS + 1), test_auprcs, label='Test AUPRC', color='green', alpha=0.7)
    axes[1].axhline(y=val_auprc, color='red', linestyle='--', label=f'Final Val AUPRC: {val_auprc:.4f}')
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('AUPRC')
    axes[1].set_title('Test AUPRC Across Epochs')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('mlp_training_metrics.png', dpi=150)
    plt.show()
    print("\nPlot saved to 'mlp_training_metrics.png'")

    # Save the model
    torch.save(model.state_dict(), 'mlp_final_model.pth')
    print("Model saved to 'mlp_final_model.pth'")


if __name__ == "__main__":
    main()
