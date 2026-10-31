#!/usr/bin/env python3
"""
MLP loss sweep across EMBEDDING TYPES for SAM MTase substrate classification.

Same loss/epoch sweep as the SAM loss test, but the variable is the EMBEDDING
rather than the loss config. Four embedding conditions are compared:

    ProtT5            embeddings_SAM_ProtT5.h5         (sequence PLM, 1024-D)
    ProstT5           embeddings_SAM_ProstT5.h5        (structure-aware PLM, 1024-D)
    ProtT5+ProstT5    concat                            (2048-D)
    ProtT5+ESM-IF     embeddings_SAM_ESM-IF.h5         (1024 + 512 structural GNN)

Fairness rules
--------------
- Every condition is evaluated on the SAME proteins: the common intersection of
  entries present in the deduplicated SAM_swiss CSV AND all three embedding H5s
  (ESM-IF drops the ~26 entries without an AFDB structure).
- The cluster-disjoint train/val split is computed ONCE on the ProtT5 embeddings
  at threshold 0.6 and reused across all embeddings and loss types, so only the
  FEATURES vary between conditions, never the split.
- The CSV is the deduplicated one (7,405 unique proteins), so no double-counting.

For each (embedding x loss_type x iteration), the MLP (input -> 512 -> 256 -> 4)
is trained to 20,000 epochs with 1,000-epoch checkpoints, under the same
PerLabelUniformLoss sweep (bce / weighted_bce / sqrt_weighted_bce / focal).

Env: smstpred (`make env`).

Outputs (./mlp_embedding_comparison_output/)
  - mlp_embedding_comparison_summary.csv      mean +/- SD per embedding/loss/epoch
  - mlp_embedding_comparison_detailed.csv     per-iteration
  - mlp_embedding_comparison_per_label.csv    per-label
  - mlp_embedding_comparison_metrics.png      metric-vs-epoch, faceted by embedding
"""

import os
import random
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import (
    matthews_corrcoef, precision_recall_curve, auc, roc_auc_score,
    f1_score, accuracy_score,
)
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import matplotlib.pyplot as plt

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)

# ============================================================================
# Configuration
# ============================================================================
DATA_DIR = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[1] / "data"))
# Input files live in DATA_DIR under their Zenodo names (see data/README.md).
DATASETS = {
    "full":     DATA_DIR / "labels_SAM_full.csv",       # deduped: 7,405 proteins
    "balanced": DATA_DIR / "labels_SAM_balanced.csv",   # deduped: 2,924 proteins
}
SOURCES = {
    "prott5":  (DATA_DIR / "embeddings_SAM_ProtT5.h5",   "flat"),        # Entry -> 1024
    "prostt5": (DATA_DIR / "embeddings_SAM_ProstT5.h5",  "mean_group"),  # Entry/mean -> 1024
    "esmif":   (DATA_DIR / "embeddings_SAM_ESM-IF.h5",   "flat"),        # Entry -> 512
}
# Each condition is a list of sources to concatenate (in order).
CONDITIONS = {
    "ProtT5":          ["prott5"],
    "ProstT5":         ["prostt5"],
    "ESM-IF":          ["esmif"],
    "ProtT5+ESM-IF":   ["prott5", "esmif"],
    "ProstT5+ESM-IF":  ["prostt5", "esmif"],
}

LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
LOSS_TYPES = ["bce", "weighted_bce", "sqrt_weighted_bce", "focal"]
THRESHOLD = 0.6
MAX_EPOCHS = 20000
# Dense early grid (every 100 up to 1k) to resolve fast learners like ProstT5,
# which already peak by ~1k; coarse 1k spacing thereafter. Nothing learns past
# ~18k, so MAX_EPOCHS stays 20000 (no tail extension).
EPOCH_CHECKPOINTS = list(range(100, 1000, 100)) + list(range(1000, 21000, 1000))
LEARNING_RATE = 0.001
N_ITERATIONS = 6

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

OUT_DIR = Path(__file__).resolve().parent / "mlp_embedding_comparison_output"


# ============================================================================
# Embedding loading
# ============================================================================

def load_source(path, layout):
    """Return {Entry: np.float32 vector} for one embedding H5."""
    out = {}
    with h5py.File(path, "r") as f:
        if layout == "flat":
            for k in f.keys():
                out[k] = np.asarray(f[k][:], dtype=np.float32)
        elif layout == "mean_group":
            for k in f.keys():
                out[k] = np.asarray(f[k]["mean"][:], dtype=np.float32)
        else:
            raise ValueError(layout)
    return out


# ============================================================================
# Clustering (in-memory; intersection is ~7k so the matrix is small)
# ============================================================================

def greedy_clustering_from_matrix(sim_matrix, threshold):
    n = len(sim_matrix)
    order = np.argsort(-np.linalg.norm(sim_matrix, axis=1))
    clusters, labels, cid = {}, -1 * np.ones(n, dtype=int), 0
    for i in tqdm(order, desc="cluster", leave=False):
        if labels[i] != -1:
            continue
        assigned = False
        for rep in clusters:
            if sim_matrix[i, clusters[rep][0]] >= threshold:
                labels[i] = rep; clusters[rep].append(i); assigned = True; break
        if not assigned:
            clusters[cid] = [i]; labels[i] = cid; cid += 1
    return clusters


def split_clusters_train_val(clusters, y, train_ratio=0.9, max_attempts=1000):
    cids = list(clusters.keys())
    for attempt in range(max_attempts):
        random.shuffle(cids)
        k = int(train_ratio * len(cids))
        tr = [i for c in cids[:k] for i in clusters[c]]
        va = [i for c in cids[k:] for i in clusters[c]]
        if all(y[va][:, j].sum() > 0 for j in range(y.shape[1])):
            return tr, va
    raise ValueError("no valid cluster split")


def split_train_pool_for_epoch(n, train_ratio=0.9):
    idx = list(range(n)); random.shuffle(idx)
    k = int(train_ratio * n)
    return idx[:k], idx[k:]


# ============================================================================
# Model / loss / metrics  (identical to the SAM loss test)
# ============================================================================

class MultilabelNN(nn.Module):
    def __init__(self, input_size, output_size):
        super().__init__()
        self.fc1 = nn.Linear(input_size, 512)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, output_size)

    def forward(self, x):
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))  # raw logits


def compute_pos_weight(y):
    counts = y.sum(axis=0).clip(min=1)
    return torch.tensor((len(y) - counts) / counts, dtype=torch.float32)


class PerLabelUniformLoss(nn.Module):
    def __init__(self, loss_type, pos_weight, gamma=2.0):
        super().__init__()
        self.loss_type, self.gamma = loss_type, gamma
        self.register_buffer("pos_weight", pos_weight)
        self.register_buffer("sqrt_pos_weight", torch.sqrt(pos_weight))

    def forward(self, logits, targets):
        total = 0.0
        for i in range(logits.shape[1]):
            li, ti = logits[:, i], targets[:, i]
            if self.loss_type == "bce":
                total += nn.functional.binary_cross_entropy_with_logits(li, ti)
            elif self.loss_type == "weighted_bce":
                total += nn.functional.binary_cross_entropy_with_logits(li, ti, pos_weight=self.pos_weight[i:i+1])
            elif self.loss_type == "sqrt_weighted_bce":
                total += nn.functional.binary_cross_entropy_with_logits(li, ti, pos_weight=self.sqrt_pos_weight[i:i+1])
            elif self.loss_type == "focal":
                bce = nn.functional.binary_cross_entropy_with_logits(li, ti, reduction="none")
                p = torch.sigmoid(li)
                p_t = ti * p + (1 - ti) * (1 - p)
                pw = ti * self.pos_weight[i] + (1 - ti)
                total += (((1 - p_t) ** self.gamma) * bce * pw).mean()
        return total


def _macro(fn, y_true, y_pred):
    vals = [fn(y_true[:, i], y_pred[:, i]) for i in range(y_true.shape[1])
            if len(np.unique(y_true[:, i])) > 1]
    return float(np.mean(vals)) if vals else 0.0


def compute_auprc(y, p):
    return _macro(lambda t, s: auc(*precision_recall_curve(t, s)[1::-1]), y, p)


def all_metrics(y, probs, thr=0.5):
    pred = (probs >= thr).astype(int)
    return {
        "auprc": compute_auprc(y, probs),
        "auroc": _macro(roc_auc_score, y, probs),
        "mcc": _macro(matthews_corrcoef, y, pred),
        "f1": _macro(f1_score, y, pred),
        "accuracy": float(np.mean([accuracy_score(y[:, i], pred[:, i]) for i in range(y.shape[1])])),
    }


def per_label_metrics(y, probs, labels, thr=0.5):
    pred = (probs >= thr).astype(int)
    out = {}
    for i, lab in enumerate(labels):
        if len(np.unique(y[:, i])) > 1:
            pr, rc, _ = precision_recall_curve(y[:, i], probs[:, i])
            out[lab] = {"auprc": auc(rc, pr), "auroc": roc_auc_score(y[:, i], probs[:, i]),
                        "mcc": matthews_corrcoef(y[:, i], pred[:, i]),
                        "f1": f1_score(y[:, i], pred[:, i]),
                        "accuracy": accuracy_score(y[:, i], pred[:, i])}
        else:
            out[lab] = {k: np.nan for k in ["auprc", "auroc", "mcc", "f1", "accuracy"]}
    return out


def find_mcc_optimal_thresholds(y_true, probs):
    """Per-label MCC-optimal thresholds tuned on the given (TRAIN) predictions,
    matching the threshold optimization used in the classifier comparison."""
    thr = np.full(len(LABELS), 0.5)
    for i in range(len(LABELS)):
        if len(np.unique(y_true[:, i])) < 2:
            continue
        best_mcc, best_t = -2.0, 0.5
        for t in np.unique(np.percentile(probs[:, i], np.linspace(0, 100, 200))):
            m = matthews_corrcoef(y_true[:, i], (probs[:, i] >= t).astype(int))
            if m > best_mcc:
                best_mcc, best_t = m, t
        thr[i] = best_t
    return thr


def train_mlp_with_checkpoints(X_tr, y_tr, X_val, y_val, loss_type):
    torch.set_num_threads(20)
    model = MultilabelNN(X_tr.shape[1], len(LABELS)).to(DEVICE)
    criterion = PerLabelUniformLoss(loss_type, compute_pos_weight(y_tr), gamma=2.0).to(DEVICE)
    opt = optim.Adam(model.parameters(), lr=LEARNING_RATE)

    # Pre-load full train/val tensors onto the device ONCE. The per-epoch 90%
    # resample is then an on-device randperm index -- no per-epoch numpy
    # indexing or host->device transfer (the real CPU bottleneck).
    Xtr_t = torch.tensor(X_tr, dtype=torch.float32, device=DEVICE)
    ytr_t = torch.tensor(y_tr, dtype=torch.float32, device=DEVICE)
    Xv_t = torch.tensor(X_val, dtype=torch.float32, device=DEVICE)
    n = X_tr.shape[0]
    k = max(1, int(0.9 * n))

    ckpts = {}
    cset = set(EPOCH_CHECKPOINTS)
    for epoch in tqdm(range(MAX_EPOCHS), desc=f"  MLP({loss_type})", leave=False, mininterval=2.0):
        idx = torch.randperm(n, device=DEVICE)[:k]
        model.train(); opt.zero_grad()
        criterion(model(Xtr_t[idx]), ytr_t[idx]).backward(); opt.step()
        if (epoch + 1) in cset:
            model.eval()
            with torch.no_grad():
                val_probs = torch.sigmoid(model(Xv_t)).cpu().numpy()
                train_probs = torch.sigmoid(model(Xtr_t)).cpu().numpy()
            # MCC-optimal per-label thresholds tuned on TRAIN, applied to val
            thr = find_mcc_optimal_thresholds(y_tr, train_probs)
            ckpts[epoch + 1] = {"metrics": all_metrics(y_val, val_probs, thr),
                                "per_label": per_label_metrics(y_val, val_probs, LABELS, thr)}
    return ckpts


# ============================================================================
# Main
# ============================================================================

def dump(all_rows, per_label_rows):
    """Write the running detailed + per-label CSVs (called incrementally so
    partial results survive a crash and can be inspected mid-run)."""
    pd.DataFrame(all_rows).to_csv(OUT_DIR / "mlp_embedding_comparison_detailed.csv", index=False)
    pd.DataFrame(per_label_rows).to_csv(OUT_DIR / "mlp_embedding_comparison_per_label.csv", index=False)


def run_dataset(ds_name, csv_path, src, all_rows, per_label_rows):
    """Full embedding x loss x epoch sweep for one dataset (full or balanced).
    Each dataset gets its own common intersection, clustering, and splits.
    Appends to the shared accumulators and dumps CSVs after each iteration."""
    df = pd.read_csv(csv_path, low_memory=False).drop_duplicates("Entry")
    label_map = {e: row for e, row in zip(df["Entry"], df[LABELS].values.astype(np.float32))}

    common = set(label_map)
    for name in src:
        common &= set(src[name])
    common = sorted(common)
    print(f"\n[{ds_name}] common intersection: {len(common)} proteins")

    y = np.vstack([label_map[e] for e in common]).astype(np.float32)
    feats = {name: np.vstack([src[name][e] for e in common]).astype(np.float32) for name in src}

    # Cluster ONCE on ProtT5 over this dataset's common set; reuse split across embeddings.
    print(f"[{ds_name}] clustering ProtT5 features at threshold {THRESHOLD}...")
    sim = cosine_similarity(feats["prott5"])
    clusters = greedy_clustering_from_matrix(sim, THRESHOLD)
    del sim
    print(f"[{ds_name}] clusters: {len(clusters)}")
    splits = [split_clusters_train_val(clusters, y, 0.9) for _ in range(N_ITERATIONS)]

    for cond_name, srcs in CONDITIONS.items():
        X = np.hstack([feats[s] for s in srcs]).astype(np.float32)
        print(f"\n[{ds_name}] === {cond_name}  (dim={X.shape[1]}) ===")
        for it, (tr, va) in enumerate(splits):
            for loss_type in LOSS_TYPES:
                ckpts = train_mlp_with_checkpoints(X[tr], y[tr], X[va], y[va], loss_type)
                for epoch, c in ckpts.items():
                    all_rows.append({"dataset": ds_name, "embedding": cond_name, "dim": X.shape[1],
                                     "loss_type": loss_type, "iteration": it + 1,
                                     "epochs": epoch, **c["metrics"]})
                    for lab, lm in c["per_label"].items():
                        per_label_rows.append({"dataset": ds_name, "embedding": cond_name,
                                               "loss_type": loss_type, "iteration": it + 1,
                                               "epochs": epoch, "label": lab, **lm})
            dump(all_rows, per_label_rows)   # incremental save after each iteration
            print(f"[{ds_name}] {cond_name} iter {it+1}/{N_ITERATIONS} done")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print("MLP embedding comparison (loss + epoch sweep per embedding; full & balanced)")
    print(f"device: {DEVICE}")
    print("=" * 70)

    # Load embedding sources once (shared across datasets).
    src = {name: load_source(path, layout) for name, (path, layout) in SOURCES.items()}
    for name in src:
        print(f"  {name}: {len(src[name])} embeddings")

    all_rows, per_label_rows = [], []
    for ds_name, csv_path in DATASETS.items():
        run_dataset(ds_name, csv_path, src, all_rows, per_label_rows)

    det = pd.DataFrame(all_rows)
    dump(all_rows, per_label_rows)   # final detailed + per-label CSVs
    metrics = ["mcc", "auprc", "auroc", "f1", "accuracy"]
    summ = (det.groupby(["dataset", "embedding", "loss_type", "epochs"])[metrics]
            .agg(["mean", "std"]).reset_index())
    summ.columns = ["_".join(c).strip("_") for c in summ.columns]
    summ.to_csv(OUT_DIR / "mlp_embedding_comparison_summary.csv", index=False)
    print(f"\nWrote CSVs to {OUT_DIR}")

    # Figure: rows = dataset, cols = metric; one line per embedding (best loss per epoch).
    datasets = list(DATASETS.keys())
    fig, axes = plt.subplots(len(datasets), len(metrics),
                             figsize=(4 * len(metrics), 4 * len(datasets)), squeeze=False)
    for r, ds_name in enumerate(datasets):
        for j, metric in enumerate(metrics):
            ax = axes[r][j]
            for cond_name in CONDITIONS:
                sub = det[(det.dataset == ds_name) & (det.embedding == cond_name)]
                if sub.empty:
                    continue
                g = sub.groupby(["epochs", "loss_type"])[metric].mean().reset_index()
                best = g.loc[g.groupby("epochs")[metric].idxmax()]
                ax.plot(best["epochs"], best[metric], marker="o", label=cond_name)
            ax.set_xlabel("epochs"); ax.set_ylabel(metric.upper())
            ax.set_title(f"{ds_name} / {metric.upper()}")
            ax.grid(True, alpha=0.3)
            if r == 0 and j == 0:
                ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "mlp_embedding_comparison_metrics.png", dpi=150)
    plt.close()
    print(f"Wrote {OUT_DIR / 'mlp_embedding_comparison_metrics.png'}")
    print("\n" + "=" * 70 + "\nCOMPLETE\n" + "=" * 70)


if __name__ == "__main__":
    main()
