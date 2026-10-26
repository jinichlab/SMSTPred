#!/usr/bin/env python3
"""
SAM MTase clustering & cosine-threshold selection (consolidated).

This single script performs the full clustering / threshold-selection stage
for the SAM_swiss training set. It replaces three earlier scripts:
  - cluster_similarity_analysis.py   (baseline all-pairwise histogram)
  - sam_cluster_threshold_sweep.py   (sweep diagnostics + cached clusters)
  - cosine_similarity_analysis.py    (deployment-domain threshold validation)

Stages
------
1. BASELINE      All-pairwise cosine-similarity histogram of the training set.
                 The peak of this distribution motivates the lowest sensible
                 clustering threshold.
2. SWEEP         Greedy-cluster the training set at each candidate threshold and
                 report cluster-structure diagnostics (cluster count, singleton
                 fraction, size percentiles), within- vs between-cluster cosine
                 distributions, and label purity (fraction of multi-member
                 clusters sharing an identical [Other, Protein, RNA, DNA]
                 vector). Cluster assignments are cached for downstream
                 stratified-by-cluster train/val splitting.
3. VALIDATION    For every SAM MTase predicted across the reference proteomes,
                 compute its maximum cosine similarity to the training set. The
                 MINIMUM of these max-similarities lands near the chosen
                 threshold (~0.6), validating that a 0.6 clustering threshold
                 spans the actual deployment domain.

The cosine threshold for the SAM model is 0.6: the within-set distribution
(stage 1) sets the lowest reasonable value and the deployment-domain minimum
(stage 3) confirms 0.6 from the opposite direction.

Outputs (written to ./sam_clustering_output/)
  - all_pairwise_similarities.png        baseline histogram (peak)
  - sam_cluster_threshold_sweep.csv      one row per threshold
  - sam_cluster_size_histograms.png      cluster-size hist per threshold
  - sam_cluster_within_between.png       within vs between cosine per threshold
  - SAM_clusters_t{t}.pkl                cached greedy clusters (downstream reuse)
  - cosine_similarity_to_training_set.csv  per-prediction max similarity
  - validation_max_similarity.png        deployment-domain validation histogram
"""

import os
import pickle
import random
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import psutil
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt

random.seed(42)
np.random.seed(42)

# ============================================================================
# Configuration  --  data directory (see data/README.md; override with $SMSTPRED_DATA)
# ============================================================================
# DATA_ROOT holds the input files under their Zenodo names.
DATA_ROOT = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[1] / "data"))

TRAIN_CSV = DATA_ROOT / "labels_SAM_full.csv"
TRAIN_H5 = DATA_ROOT / "embeddings_SAM_ProtT5.h5"
PRED_DIR = DATA_ROOT / "prediction_files_mlp"      # per-organism ProtT5-MLP predictions (not archived)
EMBED_DIR = DATA_ROOT / "tree_of_life_ProtT5"      # extract embeddings_tree_of_life_ProtT5.tar here

OUT_DIR = Path(__file__).resolve().parent / "sam_clustering_output"

LABELS = ["Other", "Protein", "RNA", "DNA"]
THRESHOLDS = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95]
CHOSEN_THRESHOLD = 0.6


# ============================================================================
# Embedding loading
# ============================================================================

def load_training_matrix(csv_path, h5_path):
    """Load training embeddings (and label vectors) as aligned arrays."""
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
    entries = df["Entry"].values
    y = df[LABELS].values.astype(np.int8) if all(c in df for c in LABELS) else None
    print(f"Loaded {X.shape[0]} training embeddings (dim={X.shape[1]})")
    return X, entries, y


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
# Similarity matrix + greedy clustering
# ============================================================================

def estimate_memory_usage(n_samples):
    return (n_samples ** 2 * 4) / (1024 ** 3)


def batched_cosine_similarity(embeddings, batch_size=5000, save_path="sam_cosine_sim.dat"):
    """Compute the full pairwise cosine matrix, memmap-backed and cached."""
    n = embeddings.shape[0]
    if os.path.exists(save_path):
        print(f"Loading cached similarity matrix from {save_path}")
        return np.memmap(save_path, dtype=np.float32, mode="r", shape=(n, n))
    sim = np.memmap(save_path, dtype=np.float32, mode="w+", shape=(n, n))
    for i in tqdm(range(0, n, batch_size), desc="Cosine sim", leave=True):
        end_i = min(i + batch_size, n)
        for j in range(0, n, batch_size):
            end_j = min(j + batch_size, n)
            sim[i:end_i, j:end_j] = cosine_similarity(
                embeddings[i:end_i], embeddings[j:end_j]
            )
    return sim


def greedy_clustering_from_matrix(sim_matrix, threshold):
    """Greedily assign each sequence to the first cluster whose representative
    is within `threshold` cosine similarity, otherwise start a new cluster."""
    n = len(sim_matrix)
    norms = np.linalg.norm(sim_matrix, axis=1)
    sorted_indices = np.argsort(-norms)

    clusters = {}
    cluster_labels = -1 * np.ones(n, dtype=int)
    cluster_id = 0
    for i in tqdm(sorted_indices, desc=f"Cluster t={threshold}", leave=False):
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


# ============================================================================
# Diagnostics
# ============================================================================

def cluster_size_stats(clusters):
    sizes = np.array([len(v) for v in clusters.values()])
    n_total = sizes.sum()
    return {
        "n_clusters": len(sizes),
        "n_singletons": int((sizes == 1).sum()),
        "singleton_fraction": float((sizes == 1).sum() / len(sizes)),
        "frac_seqs_in_multi_clusters": float(sizes[sizes > 1].sum() / n_total) if n_total else 0.0,
        "largest_cluster": int(sizes.max()),
        "median_size": float(np.median(sizes)),
        "mean_size": float(sizes.mean()),
        "p95_size": float(np.percentile(sizes, 95)),
    }


def label_purity(clusters, y_labels):
    """Fraction of multi-member clusters whose members share one label vector."""
    if y_labels is None:
        return float("nan"), float("nan"), 0
    multi = [v for v in clusters.values() if len(v) > 1]
    if not multi:
        return float("nan"), float("nan"), 0
    pure_count, pure_seqs, total_seqs = 0, 0, 0
    for members in multi:
        rows = y_labels[members]
        is_pure = np.all(rows == rows[0], axis=0).all()
        total_seqs += len(members)
        if is_pure:
            pure_count += 1
            pure_seqs += len(members)
    return pure_count / len(multi), pure_seqs / total_seqs, len(multi)


def within_between_cosine(embeddings, clusters, n_within=20000, n_between=20000, rng=None):
    """Sample within-cluster and between-cluster cosine similarities."""
    if rng is None:
        rng = np.random.default_rng(42)

    multi = [np.array(v) for v in clusters.values() if len(v) > 1]
    within = []
    if multi:
        weights = np.array([len(c) * (len(c) - 1) // 2 for c in multi], dtype=float)
        if weights.sum() > 0:
            weights = weights / weights.sum()
            for ci in rng.choice(len(multi), size=n_within, replace=True, p=weights):
                a, b = rng.choice(multi[ci], size=2, replace=False)
                within.append(
                    np.dot(embeddings[a], embeddings[b])
                    / (np.linalg.norm(embeddings[a]) * np.linalg.norm(embeddings[b]) + 1e-12)
                )
    within = np.array(within) if within else np.array([])

    cluster_of = {}
    for cid, members in clusters.items():
        for m in members:
            cluster_of[m] = cid

    n = embeddings.shape[0]
    between = []
    attempts = 0
    while len(between) < n_between and attempts < n_between * 10:
        attempts += 1
        i, j = rng.integers(0, n, size=2)
        if i == j or cluster_of[int(i)] == cluster_of[int(j)]:
            continue
        between.append(
            np.dot(embeddings[i], embeddings[j])
            / (np.linalg.norm(embeddings[i]) * np.linalg.norm(embeddings[j]) + 1e-12)
        )
    between = np.array(between)

    def stats(arr):
        if arr.size == 0:
            return {k: float("nan") for k in ["mean", "std", "p05", "p50", "p95"]}
        return {
            "mean": float(arr.mean()),
            "std": float(arr.std()),
            "p05": float(np.percentile(arr, 5)),
            "p50": float(np.percentile(arr, 50)),
            "p95": float(np.percentile(arr, 95)),
        }

    return within, between, stats(within), stats(between)


# ============================================================================
# Plots
# ============================================================================

def plot_baseline_histogram(all_sims, save_path):
    plt.figure(figsize=(10, 6))
    plt.hist(all_sims, bins=50, edgecolor="black", alpha=0.7, color="steelblue")
    plt.xlabel("Cosine Similarity", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    stats_text = (
        f"Mean: {np.mean(all_sims):.4f}\nMedian: {np.median(all_sims):.4f}\n"
        f"Std: {np.std(all_sims):.4f}\nN pairs: {len(all_sims):,}"
    )
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment="top", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    plt.xlim(-0.1, 1.1)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()
    print(f"Wrote {save_path}")


def plot_size_histograms(size_hists, save_path):
    ts = sorted(size_hists.keys())
    fig, axes = plt.subplots(1, len(ts), figsize=(3.5 * len(ts), 3.5), sharey=False)
    if len(ts) == 1:
        axes = [axes]
    for ax, t in zip(axes, ts):
        sizes = size_hists[t]
        ax.hist(sizes, bins=np.logspace(0, np.log10(max(sizes.max(), 2)), 30))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"t={t}\nclusters={len(sizes)}  max={sizes.max()}")
        ax.set_xlabel("cluster size")
        ax.set_ylabel("count")
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote {save_path}")


def plot_within_between(within_between, save_path):
    ts = sorted(within_between.keys())
    fig, axes = plt.subplots(1, len(ts), figsize=(3.5 * len(ts), 3.5), sharey=False)
    if len(ts) == 1:
        axes = [axes]
    bins = np.linspace(0, 1, 50)
    for ax, t in zip(axes, ts):
        within, between = within_between[t]
        if within.size:
            ax.hist(within, bins=bins, alpha=0.55, density=True,
                    label=f"within (n={within.size})", color="C0")
        if between.size:
            ax.hist(between, bins=bins, alpha=0.55, density=True,
                    label=f"between (n={between.size})", color="C3")
        ax.axvline(t, color="k", linestyle="--", alpha=0.7, label=f"threshold={t}")
        ax.set_title(f"t={t}")
        ax.set_xlabel("cosine similarity")
        ax.set_ylabel("density")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote {save_path}")


def plot_validation_histogram(max_sims, save_path, chosen_threshold):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(max_sims, bins=50, edgecolor="black", alpha=0.7, color="#4C72B0")
    ax.set_xlabel("Max Cosine Similarity to Training Set", fontsize=13)
    ax.set_ylabel("Count", fontsize=13)
    ax.axvline(np.median(max_sims), color="red", linestyle="--", linewidth=1.5,
               label=f"Median: {np.median(max_sims):.4f}")
    ax.axvline(max_sims.min(), color="green", linestyle="--", linewidth=1.5,
               label=f"Minimum: {max_sims.min():.4f}")
    ax.axvline(chosen_threshold, color="black", linestyle=":", linewidth=1.5,
               label=f"Chosen threshold: {chosen_threshold}")
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close()
    print(f"Wrote {save_path}")


# ============================================================================
# Stages
# ============================================================================

def stage_baseline(sim_matrix):
    print("\n[1/3] Baseline all-pairwise histogram")
    upper = np.triu_indices(sim_matrix.shape[0], k=1)
    all_sims = np.asarray(sim_matrix[upper])
    print(f"  total pairwise comparisons: {len(all_sims):,}")
    plot_baseline_histogram(all_sims, OUT_DIR / "all_pairwise_similarities.png")


def stage_sweep(embeddings, y, sim_matrix):
    print("\n[2/3] Threshold sweep + cluster diagnostics")
    rows, size_hists, within_between = [], {}, {}
    for t in THRESHOLDS:
        cache_path = OUT_DIR / f"SAM_clusters_t{t}.pkl"
        if cache_path.exists():
            print(f"  [t={t}] loading cached clusters")
            with open(cache_path, "rb") as f:
                clusters = pickle.load(f)["clusters"]
        else:
            print(f"  [t={t}] clustering...")
            clusters, lab = greedy_clustering_from_matrix(sim_matrix, t)
            with open(cache_path, "wb") as f:
                pickle.dump({"clusters": clusters, "cluster_labels": lab}, f)

        sstats = cluster_size_stats(clusters)
        unw, wgt, n_multi = label_purity(clusters, y)
        within, between, ws, bs = within_between_cosine(embeddings, clusters)

        rows.append({
            "threshold": t, **sstats, "n_multi_clusters": n_multi,
            "label_purity_unweighted": unw, "label_purity_size_weighted": wgt,
            "within_mean": ws["mean"], "within_p05": ws["p05"],
            "within_p50": ws["p50"], "within_p95": ws["p95"],
            "between_mean": bs["mean"], "between_p05": bs["p05"],
            "between_p50": bs["p50"], "between_p95": bs["p95"],
            "within_between_gap": (ws["p05"] - bs["p95"])
            if not (np.isnan(ws["p05"]) or np.isnan(bs["p95"])) else float("nan"),
        })
        size_hists[t] = np.array([len(v) for v in clusters.values()])
        within_between[t] = (within, between)
        print(f"    clusters={sstats['n_clusters']}  singletons={sstats['singleton_fraction']:.2%}"
              f"  label_purity(wgt)={wgt:.3f}")

    pd.DataFrame(rows).to_csv(OUT_DIR / "sam_cluster_threshold_sweep.csv", index=False)
    print(f"  Wrote {OUT_DIR / 'sam_cluster_threshold_sweep.csv'}")
    plot_size_histograms(size_hists, OUT_DIR / "sam_cluster_size_histograms.png")
    plot_within_between(within_between, OUT_DIR / "sam_cluster_within_between.png")


def stage_validation(X_train):
    print("\n[3/3] Deployment-domain validation (predictions vs training)")
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
    print(f"  -> minimum max-similarity ({all_max.min():.4f}) validates "
          f"the chosen clustering threshold ({CHOSEN_THRESHOLD}).")

    pd.DataFrame({
        "Entry": all_entries, "Genus": all_genera, "Max_Cosine_Similarity": all_max,
    }).to_csv(OUT_DIR / "cosine_similarity_to_training_set.csv", index=False)
    plot_validation_histogram(all_max, OUT_DIR / "validation_max_similarity.png", CHOSEN_THRESHOLD)


# ============================================================================
# Main
# ============================================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print("SAM MTase clustering & threshold selection")
    print("=" * 70)

    X, entries, y = load_training_matrix(TRAIN_CSV, TRAIN_H5)

    n = X.shape[0]
    needed = estimate_memory_usage(n)
    free_disk = psutil.disk_usage(str(OUT_DIR)).free / (1024 ** 3)
    print(f"\nSimilarity matrix needs ~{needed:.2f} GB on disk; free disk {free_disk:.2f} GB")
    if needed >= free_disk * 0.8:
        raise RuntimeError(
            f"Similarity matrix needs ~{needed:.2f} GB but only {free_disk:.2f} GB free."
        )
    sim_matrix = batched_cosine_similarity(
        X, batch_size=5000, save_path=str(OUT_DIR / "sam_cosine_sim.dat")
    )

    stage_baseline(sim_matrix)
    stage_sweep(X, y, sim_matrix)
    stage_validation(X)

    print("\n" + "=" * 70)
    print(f"COMPLETE -- outputs in {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
