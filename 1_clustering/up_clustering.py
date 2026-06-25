#!/usr/bin/env python3
"""
UP enzyme clustering & cosine-threshold selection (consolidated).

UP-set analog of sam_clustering.py for the UP_swiss enzyme training set (all
reviewed Swiss-Prot enzymes, not just SAM MTases). Replaces the earlier UP-side
scripts:
  - UP_within_set_cosine_similarity.py  (streaming baseline histogram)
  - cosine_similarity_histogram.py      (older UP histogram, superseded)

Stages
------
1. BASELINE   Within-set pairwise cosine histogram. UP is large (~201k enzymes;
              a full in-memory matrix would be ~151 GB), so the baseline is
              computed by STREAMING: similarities are binned on-the-fly with
              Welford running statistics rather than stored. The point where the
              distribution drops off motivates the chosen threshold (0.5 for UP).
2. VALIDATION For every SAM MTase predicted across the reference proteomes by
              the UP-trained model, compute its maximum cosine similarity to the
              UP training set, as a deployment-domain reference.

Why there is no clustering sweep here (unlike sam_clustering.py)
---------------------------------------------------------------
A SAM-style greedy clustering sweep is infeasible at UP scale: the full
similarity matrix is ~151 GB and greedy clustering of ~201k sequences is
O(n x n_clusters). It was never part of the original UP analysis, and it is not
needed: the released UP model (UP_genus_final_mlp_model.pth) is trained with a
RANDOM 90/10 split, not a cluster-disjoint split, so UP clustering diagnostics
are not load-bearing. The UP threshold (0.5) is justified by the within-set
streaming baseline alone.

Outputs (written to ./up_clustering_output/)
  - up_within_set_baseline.png           streaming baseline histogram
  - up_within_set_baseline_stats.csv     baseline running statistics
  - up_cosine_similarity_to_training_set.csv  per-prediction max similarity
  - up_validation_max_similarity.png     deployment-domain validation histogram
"""

import random
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt

random.seed(42)
np.random.seed(42)

# ============================================================================
# Configuration  --  edit DATA_ROOT to match your layout
# ============================================================================
DATA_ROOT = Path("/home/nroethler/Code")

TRAIN_CSV = DATA_ROOT / "files" / "UP_swiss_enzymes_substrate_types.csv"
TRAIN_H5 = DATA_ROOT / "files" / "per-protein.h5"
# Deployment set: the same tree-of-life SAM MTases the SAM model predicts on
# (prediction_files_mlp), scored here against the UP training set. This gives a
# like-for-like comparison with sam_clustering.py's validation (identical
# sequences, UP training as the reference instead of SAM training).
PRED_DIR = DATA_ROOT / "classifier" / "reference_files" / "prediction_files_mlp"
EMBED_DIR = DATA_ROOT / "classifier" / "reference_files" / "embeddings"

OUT_DIR = Path(__file__).resolve().parent / "up_clustering_output"

LABELS = ["Other", "Protein", "RNA", "DNA"]
CHOSEN_THRESHOLD = 0.5


# ============================================================================
# Embedding loading
# ============================================================================

def load_training_matrix(csv_path, h5_path):
    """Load training embeddings as a matrix."""
    df = pd.read_csv(csv_path, low_memory=False)
    data_map = {}
    with h5py.File(h5_path, "r") as f:
        for seq_id, embedding in f.items():
            data_map[seq_id] = np.array(embedding)

    df["Embeddings"] = df["Entry"].map(data_map)
    df = df.dropna(subset=["Embeddings"])
    valid = df["Embeddings"].apply(
        lambda x: isinstance(x, np.ndarray) and x.ndim == 1 and x.shape[0] == 1024
    )
    df = df[valid]

    X = np.vstack(df["Embeddings"].values).astype("float32")
    print(f"Loaded {X.shape[0]} training embeddings (dim={X.shape[1]})")
    return X


def load_prediction_matrix(h5_path, prediction_csv):
    """Load embeddings for the entries listed in a per-organism prediction CSV."""
    pred_df = pd.read_csv(prediction_csv)
    if "Entry" not in pred_df.columns:
        return None, None
    wanted = set(pred_df["Entry"].values)

    data_map = {}
    with h5py.File(h5_path, "r") as f:
        for seq_id, embedding in f.items():
            if seq_id in wanted:
                emb = np.array(embedding)
                if emb.ndim == 1 and emb.shape[0] == 1024:
                    data_map[seq_id] = emb
    if not data_map:
        return None, None
    entries = list(data_map.keys())
    X = np.vstack([data_map[e] for e in entries]).astype("float32")
    return X, np.array(entries)


# ============================================================================
# Plots
# ============================================================================

def plot_streaming_baseline(bin_edges, hist_counts, running_count, mean_sim, median_sim, n, save_path):
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    density = hist_counts / (running_count * bin_widths)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(bin_centers, density, width=bin_widths, edgecolor="black",
           alpha=0.7, color="#4C72B0", linewidth=0.3)
    ax.set_xlabel("Cosine Similarity", fontsize=13)
    ax.set_ylabel("Density", fontsize=13)
    ax.set_title(
        f"Within-Set Pairwise Cosine Similarity\n"
        f"UP Enzymes Training Set (n={n:,} sequences, {running_count:,} pairs)",
        fontsize=14,
    )
    ax.axvline(median_sim, color="red", linestyle="--", linewidth=1.5,
               label=f"Median (approx): {median_sim:.4f}")
    ax.axvline(mean_sim, color="orange", linestyle="--", linewidth=1.5,
               label=f"Mean: {mean_sim:.4f}")
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close()
    print(f"Wrote {save_path}")


def plot_validation_histogram(max_sims, save_path, chosen_threshold):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(max_sims, bins=50, edgecolor="black", alpha=0.7, color="#4C72B0")
    ax.set_xlabel("Max Cosine Similarity to UP Training Set", fontsize=13)
    ax.set_ylabel("Count", fontsize=13)
    ax.set_title(
        "Deployment-domain reference (UP-trained model)\n"
        f"(n={len(max_sims):,} predicted SAM MTases vs UP training set)",
        fontsize=14,
    )
    ax.axvline(np.median(max_sims), color="red", linestyle="--", linewidth=1.5,
               label=f"Median: {np.median(max_sims):.4f}")
    ax.axvline(max_sims.min(), color="green", linestyle="--", linewidth=1.5,
               label=f"Minimum: {max_sims.min():.4f}")
    ax.axvline(chosen_threshold, color="black", linestyle=":", linewidth=1.5,
               label=f"Benchmark-CV threshold: {chosen_threshold}")
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close()
    print(f"Wrote {save_path}")


# ============================================================================
# Stages
# ============================================================================

def stage_baseline(X):
    """Streaming within-set pairwise cosine histogram (avoids storing the
    full ~151 GB matrix in memory)."""
    print("\n[1/2] Streaming baseline within-set histogram")
    n = X.shape[0]
    n_bins = 200
    bin_edges = np.linspace(-1.0, 1.0, n_bins + 1)
    hist_counts = np.zeros(n_bins, dtype=np.int64)

    running_count = 0
    running_sum = 0.0
    running_sum_sq = 0.0
    running_min, running_max = np.inf, -np.inf

    batch_size = 1000
    for i in tqdm(range(0, n, batch_size), desc="Row batches"):
        end_i = min(i + batch_size, n)
        sim_block = cosine_similarity(X[i:end_i], X[i:])  # (end_i - i, n - i)
        for local_row in range(end_i - i):
            start_col = local_row + 1
            if start_col >= sim_block.shape[1]:
                continue
            row_sims = sim_block[local_row, start_col:]
            hist_counts += np.histogram(row_sims, bins=bin_edges)[0]
            running_count += len(row_sims)
            running_sum += row_sims.sum()
            running_sum_sq += (row_sims.astype(np.float64) ** 2).sum()
            running_min = min(running_min, row_sims.min())
            running_max = max(running_max, row_sims.max())

    mean_sim = running_sum / running_count
    std_sim = np.sqrt(running_sum_sq / running_count - mean_sim ** 2)
    cumulative = np.cumsum(hist_counts)
    median_bin = np.searchsorted(cumulative, running_count / 2)
    median_sim = (bin_edges[median_bin] + bin_edges[median_bin + 1]) / 2

    print(f"  pairs={running_count:,}  mean={mean_sim:.4f}  median~{median_sim:.4f}"
          f"  min={running_min:.4f}  max={running_max:.4f}")

    pd.DataFrame([{
        "n_sequences": n, "n_pairs": running_count, "mean": mean_sim, "std": std_sim,
        "median_approx": median_sim, "min": running_min, "max": running_max,
    }]).to_csv(OUT_DIR / "up_within_set_baseline_stats.csv", index=False)
    plot_streaming_baseline(bin_edges, hist_counts, running_count, mean_sim,
                            median_sim, n, OUT_DIR / "up_within_set_baseline.png")


def stage_validation(X_train):
    print("\n[2/2] Deployment-domain reference (UP predictions vs UP training)")
    if not PRED_DIR.exists():
        print(f"  prediction dir {PRED_DIR} not found -- skipping validation stage")
        return
    pred_files = sorted(PRED_DIR.glob("*_predictions.csv"))
    print(f"  processing {len(pred_files)} proteome prediction files")

    all_max, all_entries, all_genera = [], [], []
    batch_size = 500
    for pred_file in tqdm(pred_files, desc="Proteomes"):
        genus = pred_file.name.replace("_predictions.csv", "")
        h5_path = EMBED_DIR / f"{genus}_embeddings.h5"
        if not h5_path.exists():
            continue
        X_pred, pred_entries = load_prediction_matrix(h5_path, pred_file)
        if X_pred is None:
            continue
        max_sims = np.zeros(X_pred.shape[0], dtype=np.float32)
        for i in range(0, X_pred.shape[0], batch_size):
            block = cosine_similarity(X_pred[i:i + batch_size], X_train)
            max_sims[i:i + batch_size] = block.max(axis=1)
        all_max.extend(max_sims.tolist())
        all_entries.extend(pred_entries.tolist())
        all_genera.extend([genus] * len(pred_entries))

    if not all_max:
        print("  no prediction embeddings found -- skipping validation outputs")
        return
    all_max = np.array(all_max)
    print(f"  predicted sequences analyzed: {len(all_max):,}")
    print(f"  min={all_max.min():.4f}  median={np.median(all_max):.4f}  mean={all_max.mean():.4f}")

    pd.DataFrame({
        "Entry": all_entries, "Genus": all_genera, "Max_Cosine_Similarity": all_max,
    }).to_csv(OUT_DIR / "up_cosine_similarity_to_training_set.csv", index=False)
    plot_validation_histogram(all_max, OUT_DIR / "up_validation_max_similarity.png", CHOSEN_THRESHOLD)


# ============================================================================
# Main
# ============================================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print("UP enzyme clustering & threshold selection")
    print("=" * 70)

    X = load_training_matrix(TRAIN_CSV, TRAIN_H5)

    # Skip the (slow) streaming baseline if it already exists from a prior run.
    if (OUT_DIR / "up_within_set_baseline.png").exists() and \
       (OUT_DIR / "up_within_set_baseline_stats.csv").exists():
        print("\n[1/2] baseline outputs already exist -- skipping")
    else:
        stage_baseline(X)

    stage_validation(X)

    print("\n" + "=" * 70)
    print(f"COMPLETE -- outputs in {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
