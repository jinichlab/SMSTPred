"""
Analyze within-cluster and between-cluster cosine similarity distributions across clustering thresholds.

This script:
1. Loads sequence embeddings
2. Computes full pairwise cosine similarity matrix
3. Clusters at thresholds 0.2 to 1.0 (step 0.05)
4. Creates histograms of within-cluster cosine similarities for each threshold
5. Creates histograms of between-cluster (outside) cosine similarities for each threshold
6. Creates a baseline histogram of all pairwise similarities
"""

import os
import numpy as np
import pandas as pd
import h5py
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import matplotlib.pyplot as plt


# ============================================================================
# Data Loading (from original script)
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
# Clustering Functions (from original script)
# ============================================================================

def batched_cosine_similarity_inmemory(embeddings, batch_size=5000):
    """Compute cosine similarity in batches, store in memory."""
    n = embeddings.shape[0]
    sim_matrix = np.zeros((n, n), dtype=np.float32)

    for i in tqdm(range(0, n, batch_size), desc="Computing Cosine Similarity", leave=True):
        end_i = min(i + batch_size, n)
        for j in range(0, n, batch_size):
            end_j = min(j + batch_size, n)
            sim_matrix[i:end_i, j:end_j] = cosine_similarity(
                embeddings[i:end_i], embeddings[j:end_j]
            )
    return sim_matrix


def greedy_clustering_from_matrix(sim_matrix, threshold):
    """Greedy clustering using a precomputed cosine similarity matrix."""
    n = len(sim_matrix)
    embedding_norms = np.linalg.norm(sim_matrix, axis=1)
    sorted_indices = np.argsort(-embedding_norms)

    clusters = {}
    cluster_labels = -1 * np.ones(n, dtype=int)
    cluster_id = 0

    for i in tqdm(sorted_indices, desc=f"Clustering (threshold={threshold:.2f})", leave=False):
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


# ============================================================================
# Within-Cluster Similarity Extraction
# ============================================================================

def get_within_cluster_similarities(sim_matrix, clusters):
    """
    Extract all pairwise cosine similarities within each cluster.

    Returns a flat array of all within-cluster similarity values.
    """
    within_sims = []

    for cluster_id, indices in clusters.items():
        if len(indices) < 2:
            # Skip singleton clusters (no pairs)
            continue

        # Get all pairwise similarities within this cluster
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                sim = sim_matrix[indices[i], indices[j]]
                within_sims.append(sim)

    return np.array(within_sims)


def get_between_cluster_similarities(sim_matrix, clusters, cluster_labels, sample_size=None):
    """
    Extract cosine similarities between sequences in DIFFERENT clusters.

    For each sequence, gets similarities to all sequences outside its cluster.
    Only counts each pair once (i < j).

    Args:
        sim_matrix: Full pairwise similarity matrix
        clusters: Dictionary mapping cluster_id -> list of indices
        cluster_labels: Array mapping each sequence index to its cluster
        sample_size: If set, randomly sample this many pairs (for memory efficiency)

    Returns a flat array of all between-cluster similarity values.
    """
    n = sim_matrix.shape[0]
    between_sims = []

    # For efficiency, iterate through upper triangle only
    for i in tqdm(range(n), desc="Extracting between-cluster similarities", leave=False):
        cluster_i = cluster_labels[i]
        for j in range(i + 1, n):
            cluster_j = cluster_labels[j]
            # Only include if sequences are in different clusters
            if cluster_i != cluster_j:
                between_sims.append(sim_matrix[i, j])

    between_sims = np.array(between_sims)

    # Optionally sample if there are too many pairs
    if sample_size is not None and len(between_sims) > sample_size:
        indices = np.random.choice(len(between_sims), sample_size, replace=False)
        between_sims = between_sims[indices]

    return between_sims


def get_all_pairwise_similarities(sim_matrix):
    """
    Extract all pairwise cosine similarities (upper triangle, excluding diagonal).
    """
    n = sim_matrix.shape[0]
    # Get upper triangle indices (excluding diagonal)
    upper_tri_indices = np.triu_indices(n, k=1)
    return sim_matrix[upper_tri_indices]


# ============================================================================
# Histogram Plotting
# ============================================================================

def plot_histogram(similarities, title, output_path, threshold=None, color='steelblue', empty_msg='No pairs'):
    """Create and save a histogram of cosine similarities."""
    plt.figure(figsize=(10, 6))

    if len(similarities) == 0:
        plt.text(0.5, 0.5, empty_msg,
                 ha='center', va='center', fontsize=14, transform=plt.gca().transAxes)
        plt.title(title)
    else:
        plt.hist(similarities, bins=50, edgecolor='black', alpha=0.7, color=color)
        plt.xlabel('Cosine Similarity', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.title(title, fontsize=14)

        # Add statistics
        mean_sim = np.mean(similarities)
        median_sim = np.median(similarities)
        std_sim = np.std(similarities)

        stats_text = f'Mean: {mean_sim:.4f}\nMedian: {median_sim:.4f}\nStd: {std_sim:.4f}\nN pairs: {len(similarities):,}'
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
                 fontsize=10, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Add threshold line if provided
        if threshold is not None:
            plt.axvline(x=threshold, color='red', linestyle='--', linewidth=2,
                       label=f'Clustering threshold: {threshold}')
            plt.legend(loc='upper right')

        plt.xlim(-0.1, 1.1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_combined_histogram(within_sims, between_sims, title, output_path, threshold=None):
    """Create and save a combined histogram comparing within vs between cluster similarities."""
    plt.figure(figsize=(12, 7))

    has_within = len(within_sims) > 0
    has_between = len(between_sims) > 0

    if not has_within and not has_between:
        plt.text(0.5, 0.5, 'No pairs to display',
                 ha='center', va='center', fontsize=14, transform=plt.gca().transAxes)
        plt.title(title)
    else:
        # Plot both histograms overlaid
        if has_within:
            plt.hist(within_sims, bins=50, edgecolor='black', alpha=0.6, color='steelblue',
                    label=f'Within clusters (n={len(within_sims):,})', density=True)
        if has_between:
            plt.hist(between_sims, bins=50, edgecolor='black', alpha=0.6, color='coral',
                    label=f'Between clusters (n={len(between_sims):,})', density=True)

        plt.xlabel('Cosine Similarity', fontsize=12)
        plt.ylabel('Density', fontsize=12)
        plt.title(title, fontsize=14)

        # Add statistics text
        stats_lines = []
        if has_within:
            stats_lines.append(f'Within - Mean: {np.mean(within_sims):.4f}, Std: {np.std(within_sims):.4f}')
        if has_between:
            stats_lines.append(f'Between - Mean: {np.mean(between_sims):.4f}, Std: {np.std(between_sims):.4f}')

        stats_text = '\n'.join(stats_lines)
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
                 fontsize=10, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Add threshold line if provided
        if threshold is not None:
            plt.axvline(x=threshold, color='red', linestyle='--', linewidth=2,
                       label=f'Clustering threshold: {threshold}')

        plt.legend(loc='upper right')
        plt.xlim(-0.1, 1.1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


# ============================================================================
# Main
# ============================================================================

def main():
    # Create output directory
    output_dir = 'cluster_figures'
    os.makedirs(output_dir, exist_ok=True)

    # ========================================================================
    # Load Data
    # ========================================================================
    print("=" * 70)
    print("Cluster Similarity Analysis")
    print("=" * 70)

    print("\nLoading dataset...")
    SAM_swiss_df = pd.read_csv('../files/SAM_swiss_enzymes_substrate_types.csv')
    print(f"Original dataset shape: {SAM_swiss_df.shape}")

    SAM_swiss_prott5_file = '../files/SAM_swiss.h5'
    SAM_swiss_df = get_embeds(SAM_swiss_prott5_file, SAM_swiss_df)
    print(f"Dataset shape after loading embeddings: {SAM_swiss_df.shape}")

    # Extract embeddings
    embeddings = np.vstack(SAM_swiss_df['Embeddings'].values).astype(np.float32)
    print(f"Embeddings shape: {embeddings.shape}")

    # ========================================================================
    # Compute Full Similarity Matrix
    # ========================================================================
    print("\nComputing full cosine similarity matrix...")
    sim_matrix = batched_cosine_similarity_inmemory(embeddings, batch_size=5000)
    print(f"Similarity matrix shape: {sim_matrix.shape}")

    # ========================================================================
    # Create Baseline Histogram (All Pairwise Similarities)
    # ========================================================================
    print("\nCreating baseline histogram of all pairwise similarities...")
    all_sims = get_all_pairwise_similarities(sim_matrix)
    print(f"Total pairwise comparisons: {len(all_sims):,}")

    plot_histogram(
        all_sims,
        'All Pairwise Cosine Similarities (Before Clustering)',
        os.path.join(output_dir, 'all_pairwise_similarities.png')
    )
    print(f"  Saved: {output_dir}/all_pairwise_similarities.png")

    # ========================================================================
    # Cluster at Each Threshold and Create Histograms
    # ========================================================================
    thresholds = np.arange(0.2, 1.01, 0.05)  # 0.2 to 1.0 inclusive

    # Store summary statistics
    summary_data = []

    print(f"\nClustering at {len(thresholds)} thresholds and creating histograms...")

    for threshold in tqdm(thresholds, desc="Processing thresholds"):
        threshold = round(threshold, 2)  # Avoid floating point issues

        # Cluster at this threshold
        clusters, cluster_labels = greedy_clustering_from_matrix(sim_matrix, threshold)
        n_clusters = len(clusters)

        # Get within-cluster similarities
        within_sims = get_within_cluster_similarities(sim_matrix, clusters)

        # Get between-cluster similarities
        between_sims = get_between_cluster_similarities(sim_matrix, clusters, cluster_labels)

        # Calculate statistics
        n_singletons = sum(1 for c in clusters.values() if len(c) == 1)
        n_non_singletons = n_clusters - n_singletons

        # Within-cluster statistics
        if len(within_sims) > 0:
            within_mean = np.mean(within_sims)
            within_median = np.median(within_sims)
            within_std = np.std(within_sims)
            within_min = np.min(within_sims)
            within_max = np.max(within_sims)
        else:
            within_mean = within_median = within_std = within_min = within_max = np.nan

        # Between-cluster statistics
        if len(between_sims) > 0:
            between_mean = np.mean(between_sims)
            between_median = np.median(between_sims)
            between_std = np.std(between_sims)
            between_min = np.min(between_sims)
            between_max = np.max(between_sims)
        else:
            between_mean = between_median = between_std = between_min = between_max = np.nan

        summary_data.append({
            'threshold': threshold,
            'n_clusters': n_clusters,
            'n_singletons': n_singletons,
            'n_non_singleton_clusters': n_non_singletons,
            # Within-cluster stats
            'n_within_pairs': len(within_sims),
            'within_mean_similarity': within_mean,
            'within_median_similarity': within_median,
            'within_std_similarity': within_std,
            'within_min_similarity': within_min,
            'within_max_similarity': within_max,
            # Between-cluster stats
            'n_between_pairs': len(between_sims),
            'between_mean_similarity': between_mean,
            'between_median_similarity': between_median,
            'between_std_similarity': between_std,
            'between_min_similarity': between_min,
            'between_max_similarity': between_max,
        })

        # Create within-cluster histogram
        title = f'Within-Cluster Cosine Similarities (Threshold = {threshold:.2f})\n{n_clusters} clusters ({n_singletons} singletons)'
        output_path = os.path.join(output_dir, f'within_cluster_sims_threshold_{threshold:.2f}.png')
        plot_histogram(within_sims, title, output_path, threshold=threshold, color='steelblue',
                      empty_msg='No within-cluster pairs\n(all singletons)')

        # Create between-cluster histogram
        title = f'Between-Cluster Cosine Similarities (Threshold = {threshold:.2f})\n{n_clusters} clusters ({n_singletons} singletons)'
        output_path = os.path.join(output_dir, f'between_cluster_sims_threshold_{threshold:.2f}.png')
        plot_histogram(between_sims, title, output_path, threshold=threshold, color='coral',
                      empty_msg='No between-cluster pairs\n(single cluster)')

        # Create combined histogram (within vs between)
        title = f'Within vs Between Cluster Similarities (Threshold = {threshold:.2f})\n{n_clusters} clusters ({n_singletons} singletons)'
        output_path = os.path.join(output_dir, f'combined_sims_threshold_{threshold:.2f}.png')
        plot_combined_histogram(within_sims, between_sims, title, output_path, threshold=threshold)

    # ========================================================================
    # Save Summary Statistics
    # ========================================================================
    summary_df = pd.DataFrame(summary_data)
    summary_path = os.path.join(output_dir, 'clustering_summary_statistics.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSummary statistics saved to: {summary_path}")

    # Print summary table
    print("\n" + "=" * 130)
    print("Summary Statistics Across Thresholds")
    print("=" * 130)
    print(f"{'Threshold':<10} {'Clusters':<10} {'Singletons':<11} {'Within Pairs':<14} {'Within Mean':<12} {'Between Pairs':<15} {'Between Mean':<12}")
    print("-" * 130)
    for row in summary_data:
        within_mean_str = f"{row['within_mean_similarity']:.4f}" if not np.isnan(row['within_mean_similarity']) else "N/A"
        between_mean_str = f"{row['between_mean_similarity']:.4f}" if not np.isnan(row['between_mean_similarity']) else "N/A"
        print(f"{row['threshold']:<10.2f} {row['n_clusters']:<10} {row['n_singletons']:<11} "
              f"{row['n_within_pairs']:<14,} {within_mean_str:<12} {row['n_between_pairs']:<15,} {between_mean_str:<12}")

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"\nOutput saved to: {output_dir}/")
    print("  - all_pairwise_similarities.png (baseline)")
    print("  - within_cluster_sims_threshold_X.XX.png (17 files)")
    print("  - between_cluster_sims_threshold_X.XX.png (17 files)")
    print("  - combined_sims_threshold_X.XX.png (17 files)")
    print("  - clustering_summary_statistics.csv")
    print("=" * 70)


if __name__ == "__main__":
    main()
