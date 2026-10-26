#!/usr/bin/env python3
"""
ProstT5 MLP loss/epoch comparison for the ALL-UniProt (UP) enzyme set.

Companion to mlp_loss_comparison.py (which runs on ProtT5 embeddings). This is
the ProstT5 version for the UP dataset, so the SAM-vs-UP classifier comparison is
apples-to-apples on the embedding the deployment model actually uses (ProstT5).

What it does
------------
Under cluster-disjoint cross-validation (6 iterations), the MLP
(1024 -> 512 -> 256 -> 4) is trained up to 20,000 epochs and evaluated at a
DENSE-EARLY checkpoint grid (every 100 to 1k, then every 1k), for each of four
per-label-uniform losses: bce / weighted_bce / sqrt_weighted_bce / focal.
The dense-early grid resolves ProstT5's fast learning (it peaks by ~1k epochs).
The payoff is the PER-LABEL CSV: which loss wins for each label.

Convention (matches train_deployment_ensemble.py)
-------------------------------------------------
  - CLUSTER on ProtT5 @ 0.5 : reuse the cached UP greedy-cosine clustering
    (clusters_UP_ProtT5_t0.5.csv: 201,339 Entry -> cluster rows, aligned to
    labels_UP.csv, so the cluster-disjoint splits match the rest of the project).
  - TRAIN + EVAL on ProstT5 : features from embeddings_UP_ProstT5.h5 (/Entry/mean).
Proteins without a ProstT5 embedding are dropped (their cluster id is still known
from the cluster table, but they carry no features).

Outputs (./UP_prostt5_loss_comparison_output/)
  - UP_prostt5_mlp_loss_comparison_summary.csv      mean +/- SD per loss/epoch
  - UP_prostt5_mlp_loss_comparison_detailed.csv     per-iteration metrics
  - UP_prostt5_mlp_loss_comparison_per_label.csv    per-label breakdown
  - UP_prostt5_mlp_loss_comparison_metrics.png
  - best_UP_prostt5_mlp_loss_comparison_model.pth

Env: classifier_env (GPU strongly recommended -- 201k proteins x 20k epochs).
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
from sklearn.metrics import (matthews_corrcoef, precision_recall_curve, auc,
                             roc_auc_score, f1_score, accuracy_score)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)

# ============================================================================
# Configuration
# ============================================================================
FILES_DIR = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[1] / "data"))
CSV = FILES_DIR / "labels_UP.csv"
PROSTT5_H5 = FILES_DIR / "embeddings_UP_ProstT5.h5"                  # /Entry/mean (1024)
CLUSTER_CSV = FILES_DIR / "clusters_UP_ProtT5_t0.5.csv"     # Entry, cluster (row-aligned to labels_UP.csv)

OUT_DIR = Path(__file__).resolve().parent / "UP_prostt5_loss_comparison_output"

LABELS = ["Other", "Protein", "RNA", "DNA"]
LOSS_TYPES = ["bce", "weighted_bce", "sqrt_weighted_bce", "focal"]
MAX_EPOCHS = 20000
# Dense early grid (ProstT5 peaks by ~1k), coarse 1k spacing thereafter.
EPOCH_CHECKPOINTS = list(range(100, 1000, 100)) + list(range(1000, 21000, 1000))
LEARNING_RATE = 0.001
N_ITERATIONS = 6
CLUSTER_THRESHOLD = 0.5     # documented; clustering is reused, not recomputed

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# ============================================================================
# Data loading: ProstT5 features + reused ProtT5 cluster assignments
# ============================================================================

def load_prostt5_mean(path):
    """Entry -> mean ProstT5 vector (1024,)."""
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            item = f[k]
            out[k] = np.asarray(item["mean"][:] if isinstance(item, h5py.Group) else item[:],
                                dtype=np.float32)
    return out


def load_data():
    """Return X (ProstT5), y, and per-row cluster id, for UP enzymes that have a
    ProstT5 embedding. Cluster ids come from the cached ProtT5 clustering table
    (clusters_UP_ProtT5_t0.5.csv), whose rows are aligned to the label table."""
    df = pd.read_csv(CSV, low_memory=False).drop_duplicates("Entry").reset_index(drop=True)
    cl = pd.read_csv(CLUSTER_CSV)
    assert list(cl["Entry"]) == list(df["Entry"]), (
        "cluster table and label table are not row-aligned (different entries/order)")
    entry_to_cluster = dict(zip(cl["Entry"].values, cl["cluster"].values))

    feats = load_prostt5_mean(PROSTT5_H5)
    keep = df["Entry"].apply(lambda e: e in feats)
    df = df[keep].reset_index(drop=True)
    X = np.vstack([feats[e] for e in df["Entry"]]).astype(np.float32)
    y = df[LABELS].values.astype(np.float32)
    clusters_per_row = np.array([entry_to_cluster[e] for e in df["Entry"]], dtype=int)
    return df["Entry"].values, X, y, clusters_per_row


def build_clusters_dict(clusters_per_row):
    """cluster id -> list of feature-array row indices (after dropping no-embed rows)."""
    clusters = {}
    for i, c in enumerate(clusters_per_row):
        clusters.setdefault(int(c), []).append(i)
    return clusters


# ============================================================================
# Cluster-disjoint splitting
# ============================================================================

def split_clusters_train_val(clusters, y_labels, train_ratio=0.9, max_attempts=1000):
    """Whole clusters -> train (90%) / val (10%), ensuring every label has a
    positive in validation."""
    cluster_ids = list(clusters.keys())
    n_labels = y_labels.shape[1]
    for attempt in range(max_attempts):
        random.shuffle(cluster_ids)
        train_size = int(train_ratio * len(cluster_ids))
        train_idx = [i for c in cluster_ids[:train_size] for i in clusters[c]]
        val_idx = [i for c in cluster_ids[train_size:] for i in clusters[c]]
        y_val = y_labels[val_idx]
        if all(y_val[:, j].sum() > 0 for j in range(n_labels)):
            if attempt > 0:
                print(f"    valid split after {attempt + 1} attempts")
            return train_idx, val_idx
    raise ValueError(f"no valid cluster split after {max_attempts} attempts")


# ============================================================================
# Model, loss, metrics
# ============================================================================

class MultilabelNN(nn.Module):
    def __init__(self, in_dim, out_dim):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, 512)
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, out_dim)
        self.relu = nn.ReLU()

    def forward(self, x):                  # returns logits
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))


def compute_pos_weight(y):
    counts = y.sum(axis=0).clip(min=1)
    return torch.tensor((len(y) - counts) / counts, dtype=torch.float32)


class PerLabelUniformLoss(nn.Module):
    """Same loss for each label, summed (no averaging across labels)."""
    def __init__(self, loss_type, pos_weight, gamma=2.0):
        super().__init__()
        self.loss_type = loss_type
        self.gamma = gamma
        self.register_buffer("pos_weight", pos_weight)
        self.register_buffer("sqrt_pos_weight", torch.sqrt(pos_weight))

    def forward(self, logits, targets):
        total = 0.0
        for i in range(logits.shape[1]):
            li, ti = logits[:, i], targets[:, i]
            if self.loss_type == "bce":
                total += nn.functional.binary_cross_entropy_with_logits(li, ti)
            elif self.loss_type == "weighted_bce":
                total += nn.functional.binary_cross_entropy_with_logits(
                    li, ti, pos_weight=self.pos_weight[i:i + 1])
            elif self.loss_type == "sqrt_weighted_bce":
                total += nn.functional.binary_cross_entropy_with_logits(
                    li, ti, pos_weight=self.sqrt_pos_weight[i:i + 1])
            elif self.loss_type == "focal":
                bce = nn.functional.binary_cross_entropy_with_logits(li, ti, reduction="none")
                prob = torch.sigmoid(li)
                p_t = ti * prob + (1 - ti) * (1 - prob)
                pw = ti * self.pos_weight[i] + (1 - ti)
                total += (((1 - p_t) ** self.gamma) * bce * pw).mean()
        return total


def per_label_and_overall(y_true, probs, labels, threshold=0.5):
    pred = (probs >= threshold).astype(int)
    per = {}
    for i, lab in enumerate(labels):
        if len(np.unique(y_true[:, i])) > 1:
            prec, rec, _ = precision_recall_curve(y_true[:, i], probs[:, i])
            per[lab] = {"auprc": auc(rec, prec),
                        "auroc": roc_auc_score(y_true[:, i], probs[:, i]),
                        "mcc": matthews_corrcoef(y_true[:, i], pred[:, i]),
                        "f1": f1_score(y_true[:, i], pred[:, i], zero_division=0),
                        "accuracy": accuracy_score(y_true[:, i], pred[:, i])}
        else:
            per[lab] = {k: np.nan for k in ["auprc", "auroc", "mcc", "f1", "accuracy"]}
    metr = {m: float(np.nanmean([per[l][m] for l in labels])) for m in
            ["auprc", "auroc", "mcc", "f1", "accuracy"]}
    return metr, per


# ============================================================================
# GPU training with checkpoint evaluation
# ============================================================================

def train_with_checkpoints(Xtr, ytr, Xval, yval, iteration, loss_type):
    """Xtr/ytr: train-pool (device tensors). Xval/yval: held-out (val device tensor / np)."""
    print(f"\n  Iter {iteration + 1}  loss={loss_type}  (train pool {Xtr.shape[0]}, "
          f"val {Xval.shape[0]})")
    model = MultilabelNN(Xtr.shape[1], len(LABELS)).to(DEVICE)
    pos_weight = compute_pos_weight(ytr.cpu().numpy()).to(DEVICE)
    criterion = PerLabelUniformLoss(loss_type, pos_weight).to(DEVICE)
    opt = optim.Adam(model.parameters(), lr=LEARNING_RATE)

    n = Xtr.shape[0]
    k = int(0.9 * n)
    ckpt_set = set(EPOCH_CHECKPOINTS)
    results, best = {}, {"auprc": -1, "state": None, "epoch": -1}
    train_losses, test_losses = [], []

    for epoch in range(1, MAX_EPOCHS + 1):
        perm = torch.randperm(n, device=DEVICE)
        tr, te = perm[:k], perm[k:]
        model.train()
        opt.zero_grad()
        loss = criterion(model(Xtr[tr]), ytr[tr])
        loss.backward()
        opt.step()
        train_losses.append(float(loss.detach()))

        model.eval()
        with torch.no_grad():
            lt = float(criterion(model(Xtr[te]), ytr[te]))
        test_losses.append(lt)

        if epoch in ckpt_set:
            model.eval()
            with torch.no_grad():
                probs = torch.sigmoid(model(Xval)).cpu().numpy()
            metr, per = per_label_and_overall(yval, probs, LABELS)
            results[epoch] = {"metrics": metr, "per_label": per}
            if metr["auprc"] > best["auprc"]:
                best = {"auprc": metr["auprc"],
                        "state": {kk: v.detach().cpu().clone() for kk, v in model.state_dict().items()},
                        "epoch": epoch, "loss_type": loss_type, "iteration": iteration + 1}
            print(f"    epoch {epoch:5d}: AUPRC={metr['auprc']:.4f} AUROC={metr['auroc']:.4f} "
                  f"MCC={metr['mcc']:.4f} F1={metr['f1']:.4f}")
    return results, best, train_losses, test_losses


# ============================================================================
# Main
# ============================================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print(f"UP ProstT5 MLP loss comparison  (device={DEVICE})")
    print(f"Checkpoints: {EPOCH_CHECKPOINTS}")
    print("=" * 70)

    entries, X, y, clusters_per_row = load_data()
    clusters = build_clusters_dict(clusters_per_row)
    print(f"Proteins with ProstT5 features: {len(entries)} | dim {X.shape[1]} | "
          f"{len(clusters)} clusters | label totals "
          f"{dict(zip(LABELS, y.sum(0).astype(int)))}")

    Xfull = torch.tensor(X, dtype=torch.float32, device=DEVICE)
    yfull = torch.tensor(y, dtype=torch.float32, device=DEVICE)

    all_rows, per_label_rows = [], []
    overall_best = {"auprc": -1, "state": None}
    curves = {}

    for it in range(N_ITERATIONS):
        print("\n" + "=" * 70 + f"\n  ITERATION {it + 1}/{N_ITERATIONS}\n" + "=" * 70)
        tr_idx, val_idx = split_clusters_train_val(clusters, y, train_ratio=0.9)
        tr_t = torch.tensor(tr_idx, device=DEVICE)
        Xtr, ytr = Xfull[tr_t], yfull[tr_t]
        Xval = Xfull[torch.tensor(val_idx, device=DEVICE)]
        yval = y[val_idx]
        print(f"  train pool {len(tr_idx)} | val {len(val_idx)} (cluster-disjoint)")

        for loss_type in LOSS_TYPES:
            results, best, tl, te_l = train_with_checkpoints(Xtr, ytr, Xval, yval, it, loss_type)
            curves[(it, loss_type)] = (tl, te_l)
            for epoch, ck in results.items():
                all_rows.append({"epochs": epoch, "iteration": it + 1,
                                 "loss_type": loss_type, **ck["metrics"]})
                for lab, lm in ck["per_label"].items():
                    per_label_rows.append({"epochs": epoch, "iteration": it + 1,
                                           "loss_type": loss_type, "label": lab, **lm})
            if best["auprc"] > overall_best["auprc"]:
                overall_best = best

    # ---- CSVs ----
    res = pd.DataFrame(all_rows)
    metrics_names = ["auprc", "auroc", "mcc", "f1", "accuracy"]
    summary = []
    for lt in LOSS_TYPES:
        sub_l = res[res["loss_type"] == lt]
        for ep in EPOCH_CHECKPOINTS:
            sub = sub_l[sub_l["epochs"] == ep]
            for m in metrics_names:
                summary.append({"LossType": lt, "Epochs": ep, "Metric": m.upper(),
                                "Mean": float(sub[m].mean()), "Std": float(sub[m].std())})
    pd.DataFrame(summary).to_csv(OUT_DIR / "UP_prostt5_mlp_loss_comparison_summary.csv", index=False)
    res.to_csv(OUT_DIR / "UP_prostt5_mlp_loss_comparison_detailed.csv", index=False)
    pd.DataFrame(per_label_rows).to_csv(
        OUT_DIR / "UP_prostt5_mlp_loss_comparison_per_label.csv", index=False)
    print(f"\nCSVs -> {OUT_DIR}")

    # ---- Metric-comparison figure (the key one) ----
    colors = {"bce": "forestgreen", "weighted_bce": "darkorange",
              "sqrt_weighted_bce": "mediumorchid", "focal": "steelblue"}
    fig, axes = plt.subplots(1, len(metrics_names), figsize=(4 * len(metrics_names), 5))
    for idx, m in enumerate(metrics_names):
        ax = axes[idx]
        for lt in LOSS_TYPES:
            sub_l = res[res["loss_type"] == lt]
            means = [sub_l[sub_l["epochs"] == ep][m].mean() for ep in EPOCH_CHECKPOINTS]
            stds = [sub_l[sub_l["epochs"] == ep][m].std() for ep in EPOCH_CHECKPOINTS]
            ax.errorbar(EPOCH_CHECKPOINTS, means, yerr=stds, marker="o", capsize=4,
                        label=lt.upper(), color=colors[lt], markersize=4)
        ax.set_xscale("log")
        ax.set_xlabel("Epochs (log)"); ax.set_ylabel(m.upper()); ax.set_title(m.upper())
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "UP_prostt5_mlp_loss_comparison_metrics.png", dpi=150)
    plt.close()

    torch.save(overall_best["state"], OUT_DIR / "best_UP_prostt5_mlp_loss_comparison_model.pth")
    print(f"\nBest: iter {overall_best.get('iteration')}, epoch {overall_best.get('epoch')}, "
          f"loss={overall_best.get('loss_type')} (AUPRC {overall_best['auprc']:.4f})")
    print("=" * 70 + "\nUP ProstT5 COMPLETE\n" + "=" * 70)


if __name__ == "__main__":
    main()
