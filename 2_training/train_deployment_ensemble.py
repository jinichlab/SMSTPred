#!/usr/bin/env python3
"""
Stage 3 / Script 1 - build the deployment ProstT5 MLP classifier + thresholds.

Produces the two things needed to label new proteins (Script 2 consumes them):

  1. Frozen per-label decision thresholds, derived from cluster-disjoint
     OUT-OF-FOLD (OOF) predictions. Each protein is scored by a model that never
     saw it (or its homologs), so the thresholds are calibrated to the
     remote-homology regime deployment actually faces - not to in-sample
     over-confidence.
  2. A 6-member deep ensemble trained on ALL 7,405 labeled proteins (the
     production model). Averaging members smooths the per-seed variance and the
     member spread doubles as a confidence signal.

Pipeline
--------
  A. Load ProstT5 embeddings (deployment features, 1024-D) + ProtT5 embeddings
     (used ONLY for clustering, to match the rest of the project) + labels, for
     the deduplicated SAM_swiss set.
  B. Greedy cosine clustering on ProtT5 @ 0.6, then pack the clusters into K=5
     folds (whole clusters per fold; balanced by protein count).
  C. OOF pass: for each fold, train an MLP on the other 4 folds and predict the
     held-out fold -> every protein gets one "as-if-new" prediction. Pool them.
  D. Per-label MCC-optimal thresholds on the pooled OOF predictions. Freeze.
  E. Train the 6-seed ensemble on ALL proteins; save every member's weights.

Outputs (./deployment_model/)
  - thresholds.json            frozen per-label cutoffs + pooled-OOF metrics
  - oof_predictions.csv        per-protein OOF scores, truth, cluster, fold
  - ensemble_member_{0..5}.pt  final model weights (trained on all data)
  - model_config.json          architecture, label order, provenance

Env: classifier_env (GPU). Config matches the CV-validated ProstT5 winner:
MLP 1024->512->256->4, 300 epochs, per-label loss {Other: sqrt_weighted_bce;
Protein/RNA/DNA: bce}.
"""

import os
import json
import random
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import matthews_corrcoef, f1_score
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)

# ============================================================================
# Configuration
# ============================================================================
DATA_ROOT = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[1] / "data"))
FILES_DIR = DATA_ROOT
FULL_CSV = FILES_DIR / "labels_SAM_full.csv"
PROSTT5_H5 = FILES_DIR / "embeddings_SAM_ProstT5.h5"     # deployment features (mean group)
PROTT5_H5 = FILES_DIR / "embeddings_SAM_ProtT5.h5"              # clustering only (flat)

OUT_DIR = Path(__file__).resolve().parent / "deployment_model"

LABELS = ["Other", "Protein", "RNA", "DNA"]         # output index order
CLUSTER_THRESHOLD = 0.6
K_FOLDS = 5
N_ENSEMBLE = 6
EPOCHS = 300
LEARNING_RATE = 0.001

# ProstT5 full-set winner: per-label mixed loss (label idx [Other,Protein,RNA,DNA]).
_SQ, _BCE, _WB, _FOC = "sqrt_weighted_bce", "bce", "weighted_bce", "focal"
LOSS_MAP = {0: _SQ, 1: _BCE, 2: _BCE, 3: _BCE}

HIDDEN = (512, 256)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# ============================================================================
# Data loading
# ============================================================================

def load_flat(path):
    with h5py.File(path, "r") as f:
        return {k: np.array(f[k], dtype=np.float32) for k in f.keys()}


def load_mean_group(path):
    with h5py.File(path, "r") as f:
        return {k: np.array(f[k]["mean"], dtype=np.float32) for k in f.keys()}


def load_data():
    """Return entries, ProstT5 feature matrix, ProtT5 (clustering) matrix, labels."""
    prostt5 = load_mean_group(PROSTT5_H5)
    prott5 = load_flat(PROTT5_H5)
    df = pd.read_csv(FULL_CSV, low_memory=False).drop_duplicates("Entry")
    keep = df["Entry"].apply(lambda e: e in prostt5 and e in prott5)
    df = df[keep].reset_index(drop=True)
    entries = df["Entry"].values
    X_prostt5 = np.vstack([prostt5[e] for e in entries]).astype(np.float32)
    X_prott5 = np.vstack([prott5[e] for e in entries]).astype(np.float32)
    y = df[LABELS].values.astype(np.float32)
    return entries, X_prostt5, X_prott5, y


# ============================================================================
# Clustering (greedy cosine @ threshold) + fold packing
# ============================================================================

def greedy_clustering(sim, threshold):
    """Seed-based greedy clustering: a protein joins the first existing cluster
    whose seed it is >= threshold similar to, else it starts a new cluster."""
    n = len(sim)
    order = np.argsort(-np.linalg.norm(sim, axis=1))
    clusters, labels, cid = {}, -np.ones(n, dtype=int), 0
    for i in tqdm(order, desc="clustering", leave=False):
        if labels[i] != -1:
            continue
        assigned = False
        for rep in clusters:
            if sim[i, clusters[rep][0]] >= threshold:
                labels[i] = rep
                clusters[rep].append(i)
                assigned = True
                break
        if not assigned:
            clusters[cid] = [i]
            labels[i] = cid
            cid += 1
    return clusters


def pack_folds(clusters, n_folds, n):
    """Assign whole clusters to folds, greedily balancing by protein count."""
    fold_members = [[] for _ in range(n_folds)]
    sizes = [0] * n_folds
    for members in sorted(clusters.values(), key=len, reverse=True):
        f = int(np.argmin(sizes))
        fold_members[f].extend(members)
        sizes[f] += len(members)
    fold_id = np.empty(n, dtype=int)
    for f, mem in enumerate(fold_members):
        for i in mem:
            fold_id[i] = f
    return fold_id


# ============================================================================
# MLP + per-label mixed loss
# ============================================================================

class MultilabelNN(nn.Module):
    def __init__(self, in_dim, out_dim, hidden=HIDDEN):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hidden[0])
        self.fc2 = nn.Linear(hidden[0], hidden[1])
        self.fc3 = nn.Linear(hidden[1], out_dim)
        self.relu = nn.ReLU()

    def forward(self, x):                 # returns logits
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))


def compute_pos_weight(y):
    counts = y.sum(axis=0).clip(min=1)
    return torch.tensor((len(y) - counts) / counts, dtype=torch.float32)


class PerLabelMixedLoss(nn.Module):
    """Per-label loss from a {label_idx: loss_type} map, summed over labels."""
    def __init__(self, pos_weight, label_loss_map, gamma=2.0):
        super().__init__()
        self.pos_weight = pos_weight
        self.label_loss_map = label_loss_map
        self.gamma = gamma

    def forward(self, logits, targets):
        total = 0.0
        for i in range(logits.shape[1]):
            lt = self.label_loss_map.get(i, "weighted_bce")
            li, ti = logits[:, i], targets[:, i]
            if lt == "focal":
                bce = nn.functional.binary_cross_entropy_with_logits(li, ti, reduction="none")
                prob = torch.sigmoid(li)
                p_t = ti * prob + (1 - ti) * (1 - prob)
                w = ti * self.pos_weight[i] + (1 - ti)
                total += (((1 - p_t) ** self.gamma) * bce * w).mean()
            elif lt == "sqrt_weighted_bce":
                total += nn.functional.binary_cross_entropy_with_logits(
                    li, ti, pos_weight=torch.sqrt(self.pos_weight[i:i + 1]))
            elif lt == "bce":
                total += nn.functional.binary_cross_entropy_with_logits(li, ti)
            else:  # weighted_bce
                total += nn.functional.binary_cross_entropy_with_logits(
                    li, ti, pos_weight=self.pos_weight[i:i + 1])
        return total


def train_mlp(X, y, loss_map, epochs, lr, seed):
    """Train one MLP on (X, y) for `epochs` full-batch-90% steps; return model."""
    torch.manual_seed(seed)
    torch.set_num_threads(8)
    model = MultilabelNN(X.shape[1], y.shape[1]).to(DEVICE)
    criterion = PerLabelMixedLoss(compute_pos_weight(y).to(DEVICE), loss_map)
    opt = optim.Adam(model.parameters(), lr=lr)
    Xt = torch.tensor(X, dtype=torch.float32, device=DEVICE)
    yt = torch.tensor(y, dtype=torch.float32, device=DEVICE)
    n = X.shape[0]
    k = max(1, int(0.9 * n))
    for _ in range(epochs):
        idx = torch.randperm(n, device=DEVICE)[:k]
        model.train()
        opt.zero_grad()
        criterion(model(Xt[idx]), yt[idx]).backward()
        opt.step()
    return model


def predict_probs(model, X):
    model.eval()
    with torch.no_grad():
        xb = torch.tensor(X, dtype=torch.float32, device=DEVICE)
        return torch.sigmoid(model(xb)).cpu().numpy()


# ============================================================================
# Thresholds + metrics
# ============================================================================

def find_optimal_thresholds(y, probs, labels):
    """Per-label MCC-optimal cutoff (200-point percentile sweep of the scores)."""
    thr = np.full(len(labels), 0.5)
    for i in range(len(labels)):
        s, t = probs[:, i], y[:, i]
        if len(np.unique(t)) < 2:
            continue
        best_m, best_t = -2.0, 0.5
        for c in np.unique(np.percentile(s, np.linspace(0, 100, 200))):
            m = matthews_corrcoef(t, (s >= c).astype(int))
            if m > best_m:
                best_m, best_t = m, c
        thr[i] = best_t
    return thr


def per_label_report(y, probs, thr, labels):
    out = {}
    for i, lab in enumerate(labels):
        pred = (probs[:, i] >= thr[i]).astype(int)
        out[lab] = {"mcc": float(matthews_corrcoef(y[:, i], pred)),
                    "f1": float(f1_score(y[:, i], pred, zero_division=0)),
                    "threshold": float(thr[i]),
                    "n_positive": int(y[:, i].sum())}
    return out


# ============================================================================
# Main
# ============================================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print(f"Deployment classifier build  (device={DEVICE})")
    print("=" * 70)

    entries, X_prostt5, X_prott5, y = load_data()
    n = len(entries)
    print(f"Proteins: {n} | ProstT5 dim {X_prostt5.shape[1]} | "
          f"label totals {dict(zip(LABELS, y.sum(0).astype(int)))}")

    # --- Applicability-domain cutoff (ProstT5 space = the feature space) ---
    # For each training protein, its nearest-neighbour cosine similarity to any
    # OTHER training protein; the 1st percentile is the AD floor. At inference a
    # new protein is "in domain" if it is at least this close to the training set.
    Sp = cosine_similarity(X_prostt5)
    np.fill_diagonal(Sp, -1.0)
    ad_cutoff = float(np.percentile(Sp.max(axis=1), 1))
    del Sp
    print(f"AD cutoff (1st-pct training NN cosine, ProstT5): {ad_cutoff:.4f}")

    # --- B. cluster on ProtT5, pack into folds ---
    print(f"\nClustering on ProtT5 @ {CLUSTER_THRESHOLD} ...")
    clusters = greedy_clustering(cosine_similarity(X_prott5), CLUSTER_THRESHOLD)
    fold_id = pack_folds(clusters, K_FOLDS, n)
    print(f"  {len(clusters)} clusters -> {K_FOLDS} folds "
          f"(fold sizes: {[int((fold_id == f).sum()) for f in range(K_FOLDS)]})")

    # --- C. OOF pass: each protein predicted by a model that never saw its fold ---
    print("\nOOF pass (cluster-disjoint):")
    oof = np.zeros((n, len(LABELS)), dtype=np.float32)
    for f in range(K_FOLDS):
        te = np.where(fold_id == f)[0]
        tr = np.where(fold_id != f)[0]
        model = train_mlp(X_prostt5[tr], y[tr], LOSS_MAP, EPOCHS, LEARNING_RATE, seed=f)
        oof[te] = predict_probs(model, X_prostt5[te])
        miss = [LABELS[i] for i in range(len(LABELS)) if y[tr][:, i].sum() == 0]
        print(f"  fold {f}: train {len(tr)}  held-out {len(te)}"
              + (f"  [WARN train missing {miss}]" if miss else ""))

    # --- D. thresholds on pooled OOF ---
    thr = find_optimal_thresholds(y, oof, LABELS)
    report = per_label_report(y, oof, thr, LABELS)
    overall_mcc = float(np.mean([report[l]["mcc"] for l in LABELS]))
    print("\nFrozen thresholds (cluster-disjoint OOF):")
    for lab in LABELS:
        r = report[lab]
        print(f"  {lab:<8} thr={r['threshold']:.4f}  MCC={r['mcc']:.3f}  "
              f"F1={r['f1']:.3f}  (n+={r['n_positive']})")
    print(f"  overall mean MCC = {overall_mcc:.3f}")

    # save OOF predictions + thresholds
    oof_df = pd.DataFrame({"Entry": entries,
                           "cluster": _cluster_ids(clusters, n),
                           "fold": fold_id})
    for i, lab in enumerate(LABELS):
        oof_df[f"score_{lab}"] = oof[:, i]
        oof_df[f"true_{lab}"] = y[:, i].astype(int)
    oof_df.to_csv(OUT_DIR / "oof_predictions.csv", index=False)

    json.dump({"labels": LABELS,
               "thresholds": {lab: float(thr[i]) for i, lab in enumerate(LABELS)},
               "oof_metrics": report,
               "oof_overall_mcc": overall_mcc,
               "clustering": {"feature": "ProtT5", "threshold": CLUSTER_THRESHOLD,
                              "n_clusters": len(clusters), "k_folds": K_FOLDS}},
              open(OUT_DIR / "thresholds.json", "w"), indent=2)
    print(f"  wrote thresholds.json + oof_predictions.csv")

    # --- E. final ensemble on ALL data ---
    print(f"\nTraining {N_ENSEMBLE}-member ensemble on all {n} proteins:")
    for s in range(N_ENSEMBLE):
        model = train_mlp(X_prostt5, y, LOSS_MAP, EPOCHS, LEARNING_RATE, seed=100 + s)
        torch.save(model.state_dict(), OUT_DIR / f"ensemble_member_{s}.pt")
        print(f"  member {s} saved")

    json.dump({"embedding": "ProstT5", "embedding_h5": PROSTT5_H5.name,
               "input_dim": int(X_prostt5.shape[1]), "hidden": list(HIDDEN),
               "labels": LABELS, "loss_map": {LABELS[i]: LOSS_MAP[i] for i in range(len(LABELS))},
               "epochs": EPOCHS, "learning_rate": LEARNING_RATE,
               "n_ensemble": N_ENSEMBLE, "trained_on": int(n),
               "ad_reference_h5": PROSTT5_H5.name, "ad_cutoff": ad_cutoff,
               "built": str(date.today())},
              open(OUT_DIR / "model_config.json", "w"), indent=2)
    print(f"\nDone. Deployment model written to {OUT_DIR}")
    print("=" * 70)


def _cluster_ids(clusters, n):
    cid = np.empty(n, dtype=int)
    for c, members in clusters.items():
        for i in members:
            cid[i] = c
    return cid


if __name__ == "__main__":
    main()
