import random
import numpy as np
import pandas as pd
import h5py
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import matplotlib.pyplot as plt

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

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


def compute_upper_triangle_similarities(embeddings, batch_size=2000):
    """
    Compute upper triangle of cosine similarity matrix in batches.
    Returns flattened array of all pairwise similarities (excluding diagonal).
    """
    n = embeddings.shape[0]
    similarities = []

    print(f"Computing pairwise cosine similarities for {n} samples...")
    print(f"Total pairs to compute: {n * (n - 1) // 2:,}")

    # Normalize embeddings for faster cosine similarity computation
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    normalized_embeddings = embeddings / norms

    # Compute upper triangle in batches
    for i in tqdm(range(0, n, batch_size), desc="Computing similarities"):
        end_i = min(i + batch_size, n)

        # For each row batch, compute similarities with all columns after it
        for j in range(i, n, batch_size):
            end_j = min(j + batch_size, n)

            # Compute similarity block
            sim_block = np.dot(normalized_embeddings[i:end_i],
                              normalized_embeddings[j:end_j].T)

            # Extract upper triangle values
            if i == j:
                # Same block - get upper triangle excluding diagonal
                rows, cols = np.triu_indices(end_i - i, k=1)
                similarities.extend(sim_block[rows, cols].flatten())
            else:
                # Different blocks - all values are valid
                similarities.extend(sim_block.flatten())

    return np.array(similarities)


def main():
    SAMPLE_SIZE = 20000
    BIN_WIDTH = 0.025

    # Load data
    print("Loading dataset...")
    UP_swiss_df = pd.read_csv('../files/UP_swiss_enzymes_substrate_types.csv')
    print(f"Original dataset shape: {UP_swiss_df.shape}")

    UP_swiss_prott5_file = '../files/per-protein.h5'
    UP_swiss_df = get_embeds(UP_swiss_prott5_file, UP_swiss_df)
    print(f"Dataset shape after loading embeddings: {UP_swiss_df.shape}")

    # Random sample
    if len(UP_swiss_df) > SAMPLE_SIZE:
        print(f"\nRandomly sampling {SAMPLE_SIZE} sequences from {len(UP_swiss_df)}...")
        sample_df = UP_swiss_df.sample(n=SAMPLE_SIZE, random_state=42)
    else:
        print(f"\nDataset has {len(UP_swiss_df)} sequences, using all of them...")
        sample_df = UP_swiss_df
        SAMPLE_SIZE = len(sample_df)

    # Extract embeddings
    embeddings = np.vstack(sample_df['Embeddings'].values).astype(np.float32)
    print(f"Embeddings shape: {embeddings.shape}")

    # Compute all pairwise cosine similarities (upper triangle only)
    similarities = compute_upper_triangle_similarities(embeddings)
    print(f"\nTotal similarity pairs computed: {len(similarities):,}")

    # Print statistics
    print(f"\nCosine Similarity Statistics:")
    print(f"  Min: {similarities.min():.4f}")
    print(f"  Max: {similarities.max():.4f}")
    print(f"  Mean: {similarities.mean():.4f}")
    print(f"  Median: {np.median(similarities):.4f}")
    print(f"  Std: {similarities.std():.4f}")

    # Create histogram
    bins = np.arange(-1.0, 1.0 + BIN_WIDTH, BIN_WIDTH)

    plt.figure(figsize=(14, 6))

    # Plot 1: Count histogram
    plt.subplot(1, 2, 1)
    counts, bin_edges, patches = plt.hist(similarities, bins=bins, edgecolor='black',
                                           linewidth=0.5, alpha=0.7, color='steelblue')
    plt.xlabel('Cosine Similarity', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.title(f'Distribution of Pairwise Cosine Similarities\n({SAMPLE_SIZE:,} sequences, {len(similarities):,} pairs)',
              fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.axvline(x=0.7, color='red', linestyle='--', linewidth=2, label='Clustering threshold (0.7)')
    plt.legend()

    # Plot 2: Density histogram (normalized)
    plt.subplot(1, 2, 2)
    plt.hist(similarities, bins=bins, edgecolor='black', linewidth=0.5,
             alpha=0.7, color='steelblue', density=True)
    plt.xlabel('Cosine Similarity', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.title(f'Density of Pairwise Cosine Similarities\n({SAMPLE_SIZE:,} sequences)', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.axvline(x=0.7, color='red', linestyle='--', linewidth=2, label='Clustering threshold (0.7)')
    plt.legend()

    plt.tight_layout()
    plt.savefig('cosine_similarity_distribution.png', dpi=150)
    plt.show()
    print("\nPlot saved to 'cosine_similarity_distribution.png'")

    # Print bin counts for reference
    print("\nBin counts (similarity range: count):")
    for i in range(len(counts)):
        if counts[i] > 0:
            print(f"  [{bin_edges[i]:.3f}, {bin_edges[i+1]:.3f}): {int(counts[i]):,}")

    # Calculate percentage above clustering threshold
    above_threshold = np.sum(similarities >= 0.7)
    pct_above = 100 * above_threshold / len(similarities)
    print(f"\nPairs with similarity >= 0.7: {above_threshold:,} ({pct_above:.2f}%)")


if __name__ == "__main__":
    main()
