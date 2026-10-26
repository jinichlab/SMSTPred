#!/usr/bin/env python3
"""
Unified SAM substrate-type classifier comparison: every method, one shared test set.

One script that, for each split mode (cluster-disjoint and random) over 16
iterations, trains and tests EVERY classifier family on the full SAM_swiss set
and scores all of them on a single shared held-out test set, so the methods are
compared like for like. An optional second arm trains on the class-balanced
SAM_swiss subset instead (BAKEOFF_DATASETS=full,balanced); it is off by default
and is not part of the manuscript.

Shared-test-set design (per split mode, per iteration)
------------------------------------------------------
  1. Split the FULL SAM set into a 90% train pool and a 10% held-out TEST set
     (cluster-disjoint or random). That test set is the common evaluation set.
  2. FULL models train on the full 90% train pool.
  3. (Optional balanced arm.) BALANCED models train on the balanced set with any
     sequence whose Entry ID appears in the test set withheld (no leakage), then
     predict on the SAME test set.
  The balanced labels are read even when that arm is off: the Profile HMM draws
  random numbers per training set, so a skipped arm still has to advance the
  random stream (replay_hmm_rng) or every later split would differ.

Classifiers
-----------
  XGBoost, Random Forest, Logistic Regression, SVM (linear / RBF / poly-3 /
  sigmoid), Profile HMM (MAFFT -> hmmbuild -> hmmsearch), MMseqs2 1-NN
  (sequence-homology baseline), Foldseek 1-NN (structural-homology baseline),
  and a family of MLPs over different embeddings (ProtT5, ProstT5, SaProt,
  ESM-IF, and the ProtT5+ESM-IF / ProstT5+ESM-IF concatenations).

  The sklearn classifiers (XGBoost / RF / LR / SVMs) are trained on ProstT5
  embeddings, matching the deployed representation, so every embedding-based
  method reads the same features. Cluster-disjoint splits are still defined by
  ProtT5 cosine similarity (see CLUSTER_SOURCE) -- the split is a property of
  the data partition, not of any model, and holding it fixed keeps all 16
  classifiers on identical train/test sequences.

  Each MLP carries a shared epoch count and a per-label mixed loss (gamma=2.0)
  chosen from the embedding sweeps (see MLP_VARIANTS): 300 epochs for every
  full-set MLP, and for each embedding and label the loss with the best MCC at
  that epoch. The per-label loss lets the rare DNA head use focal where the
  weighted-BCE family overfits it at high epochs.

  Per-class decision thresholds are tuned on TRAIN predictions (MCC-optimal) for
  every classifier, so HMM/MMseqs bit-score outputs and the probabilistic models
  are all scored at a comparable operating point (never a fixed 0.5). MCC is used
  to match the deployment pipelines and the manuscript's headline metric.

Requires external tools on PATH: mmseqs, mafft, hmmbuild, hmmsearch.

Outputs (./sam_classifier_comparison_output/, or $BAKEOFF_OUT)
  - sam_classifier_comparison_detailed.csv   one row per split/dataset/iter/classifier
  - sam_classifier_comparison_summary.csv     mean +/- SD over iterations
  - sam_classifier_comparison_per_label.csv   per-label breakdown
  - sam_classifier_comparison.png             mean-metric bar chart
"""

import os
import sys
import pickle
import random
import tempfile
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier
from sklearn.base import clone
from sklearn.metrics import (
    matthews_corrcoef, precision_recall_curve, auc, roc_auc_score,
    f1_score, accuracy_score,
)
from sklearn.metrics.pairwise import cosine_similarity
from xgboost import XGBClassifier
from joblib import Parallel, delayed
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

FULL_CSV = FILES_DIR / "labels_SAM_full.csv"
BAL_CSV = FILES_DIR / "labels_SAM_balanced.csv"

# Embedding sources keyed by short name -> (h5 path, layout, dim). "flat" = the
# H5 value is the vector itself; "mean_group" = vector lives at f[entry]["mean"].
EMB_SOURCES = {
    "prott5":  (FILES_DIR / "embeddings_SAM_ProtT5.h5",  "flat",       1024),
    "prostt5": (FILES_DIR / "embeddings_SAM_ProstT5.h5", "mean_group", 1024),
    "esmif":   (FILES_DIR / "embeddings_SAM_ESM-IF.h5",  "flat",        512),
    "saprot":  (FILES_DIR / "embeddings_SAM_SaProt.h5", "flat",       1280),
}

# Feature source for the sklearn classifiers (XGBoost/RF/LR/SVMs): ProstT5, the
# deployed representation. Separate from the clustering source below.
SKLEARN_SOURCE = "prostt5"
# Similarity space that defines the cluster-disjoint splits. Kept on ProtT5: the
# 0.6 threshold is calibrated in this space (cosine_similarity_analysis.py) and
# sam_vs_up_prostt5_comparison.py clusters the same way. Changing it would move
# every classifier's train/test partition, not just the sklearn ones.
CLUSTER_SOURCE = "prott5"

OUT_DIR = Path(__file__).resolve().parent / os.environ.get("BAKEOFF_OUT", "sam_classifier_comparison_output")

# Foldseek structural-homology baseline: {Entry}.pdb structures + the binary.
STRUCT_DIR = FILES_DIR / "sam_structures"      # AlphaFold DB models of the SAM set (not archived)
FOLDSEEK = os.environ.get("FOLDSEEK_BIN", "foldseek")        # default: foldseek on PATH

LABELS = ["Other", "Protein", "RNA", "DNA"]
SPLIT_MODES = ["cluster", "random"]
N_ITERATIONS = 16
THRESHOLD = 0.6          # clustering threshold (cluster split mode)
TRAIN_RATIO = 0.9        # 90% train pool / 10% shared test set
LEARNING_RATE = 0.001
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# MLP variants: one per embedding (singletons + structural concatenations).
# Each carries a per-dataset epoch count and a per-label mixed-loss map (label
# order [Other, Protein, RNA, DNA] = idx 0,1,2,3), read off the 6-split
# embedding x loss x epoch sweeps (mlp_embedding_comparison*.py; Fig S2) by one rule:
#   full      300 epochs for every embedding; per label, the loss with the highest mean MCC over the six
#             sweep splits at epoch 300 (fixed 0.5 cutoff). An exact tie goes to the first of bce,
#             weighted_bce, sqrt_weighted_bce, focal (one tie in the 48 cells: DNA of ProtT5+ESM-IF).
#             Exception: MLP_ProstT5 keeps the DEPLOYED map (Other = sqrt_weighted_bce), so the benchmark scores
#             the model that ships; the rule alone would pick weighted_bce for Other, a 0.001-MCC tie.
#   balanced  (optional arm, not in the manuscript) 100 epochs, same rule on the balanced-set sweep.
_SQ, _BCE, _WB, _FOC = "sqrt_weighted_bce", "bce", "weighted_bce", "focal"
MLP_VARIANTS = {
    "MLP_ProtT5": {"sources": ["prott5"],
        "full":     {"epochs": 300, "loss": {0: _BCE, 1: _WB,  2: _SQ,  3: _SQ}},
        "balanced": {"epochs": 100, "loss": {0: _WB,  1: _WB,  2: _FOC, 3: _SQ}}},
    "MLP_ProstT5": {"sources": ["prostt5"],
        "full":     {"epochs": 300, "loss": {0: _SQ,  1: _BCE, 2: _BCE, 3: _BCE}},     # = deployed model
        "balanced": {"epochs": 100, "loss": {0: _SQ,  1: _FOC, 2: _FOC, 3: _WB}}},
    "MLP_ESM-IF": {"sources": ["esmif"],
        "full":     {"epochs": 300, "loss": {0: _WB,  1: _FOC, 2: _WB,  3: _FOC}},
        "balanced": {"epochs": 100, "loss": {0: _FOC, 1: _FOC, 2: _FOC, 3: _WB}}},
    "MLP_ProtT5+ESM-IF": {"sources": ["prott5", "esmif"],
        "full":     {"epochs": 300, "loss": {0: _SQ,  1: _WB,  2: _WB,  3: _SQ}},
        "balanced": {"epochs": 100, "loss": {0: _BCE, 1: _WB,  2: _FOC, 3: _SQ}}},
    "MLP_ProstT5+ESM-IF": {"sources": ["prostt5", "esmif"],
        "full":     {"epochs": 300, "loss": {0: _SQ,  1: _WB,  2: _SQ,  3: _SQ}},
        "balanced": {"epochs": 100, "loss": {0: _WB,  1: _FOC, 2: _FOC, 3: _SQ}}},
    "MLP_SaProt": {"sources": ["saprot"],
        "full":     {"epochs": 300, "loss": {0: _SQ,  1: _SQ,  2: _FOC, 3: _BCE}},
        "balanced": {"epochs": 100, "loss": {0: _SQ,  1: _WB,  2: _FOC, 3: _SQ}}},
}

# ---------------------------------------------------------------------------
# Run options (environment variables; with none set the script runs every method on the full set):
#   BAKEOFF_DATASETS=full,balanced   also train the optional balanced arm (default: full only)
#   BAKEOFF_ONLY=mlp                 train only the MLP variants; Profile-HMM's random draws are replayed so the
#                                    train/test splits stay identical to a full run
#   BAKEOFF_MLP_VARIANTS=a,b         with BAKEOFF_ONLY=mlp: only the named variants (e.g. MLP_ProstT5)
#   BAKEOFF_OUT=<dir>                output folder (relative to this script, or absolute)
# ---------------------------------------------------------------------------
DATASETS = [d.strip() for d in os.environ.get("BAKEOFF_DATASETS", "full").split(",")]
assert DATASETS and all(d in ("full", "balanced") for d in DATASETS), f"BAKEOFF_DATASETS must list full and/or balanced, got {DATASETS}"
ONLY_MLP = os.environ.get("BAKEOFF_ONLY") == "mlp"
if os.environ.get("BAKEOFF_MLP_VARIANTS"):
    _keep = [n.strip() for n in os.environ["BAKEOFF_MLP_VARIANTS"].split(",")]
    assert all(n in MLP_VARIANTS for n in _keep), f"unknown MLP variant in {_keep}; choose from {list(MLP_VARIANTS)}"
    MLP_VARIANTS = {n: MLP_VARIANTS[n] for n in _keep}

SKLEARN_CLASSIFIERS = {
    "XGBoost": XGBClassifier(eval_metric="logloss", random_state=42, n_jobs=-1),
    "Random Forest": RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1),
    "Logistic Regression": LogisticRegression(max_iter=1000, random_state=42, n_jobs=-1),
    "SVM_Linear": SVC(kernel="linear", probability=True, random_state=42),
    "SVM_RBF": SVC(kernel="rbf", probability=True, random_state=42),
    "SVM_Polynomial": SVC(kernel="poly", degree=3, probability=True, random_state=42),
    "SVM_Sigmoidal": SVC(kernel="sigmoid", probability=True, random_state=42),
}


# ============================================================================
# Tool check
# ============================================================================

def check_tools():
    missing = [t for t in ["mmseqs", "mafft", "hmmbuild", "hmmsearch"]
               if subprocess.run(["which", t], capture_output=True).returncode != 0]
    if not (shutil.which(FOLDSEEK) or os.path.exists(FOLDSEEK)):
        missing.append(f"foldseek (expected at {FOLDSEEK}; set FOLDSEEK_BIN)")
    if not STRUCT_DIR.is_dir():
        missing.append(f"structure dir {STRUCT_DIR} (run esm_if/fetch_structures.py)")
    if missing:
        sys.exit(f"Missing required tools/inputs: {missing}\n"
                 f"Install with: conda install -c bioconda mmseqs2 mafft hmmer")


# ============================================================================
# Data loading
# ============================================================================

def load_source_map(path, layout):
    """Entry -> embedding vector for one H5 source."""
    m = {}
    with h5py.File(path, "r") as f:
        if layout == "mean_group":
            for k in f.keys():
                m[k] = np.array(f[k]["mean"], dtype=np.float32)
        else:
            for k in f.keys():
                m[k] = np.array(f[k], dtype=np.float32)
    return m


def attach_sources(df, source_maps):
    """Restrict df to entries present in EVERY source; return (df, {src: matrix})
    with each matrix row-aligned to the returned df."""
    keep = df["Entry"].apply(lambda e: all(e in sm for sm in source_maps.values()))
    df = df[keep].reset_index(drop=True)
    mats = {s: np.vstack([sm[e] for e in df["Entry"]]).astype(np.float32)
            for s, sm in source_maps.items()}
    return df, mats


# ============================================================================
# Clustering (for cluster split mode)
# ============================================================================

def estimate_memory_usage(n):
    return (n ** 2 * 4) / (1024 ** 3)


def batched_cosine_similarity(embeddings, batch_size=5000, save_path="cosine_sim_matrix.dat"):
    n = embeddings.shape[0]
    if os.path.exists(save_path):
        return np.memmap(save_path, dtype=np.float32, mode="r", shape=(n, n))
    sim = np.memmap(save_path, dtype=np.float32, mode="w+", shape=(n, n))
    for i in range(0, n, batch_size):
        end_i = min(i + batch_size, n)
        for j in range(0, n, batch_size):
            end_j = min(j + batch_size, n)
            sim[i:end_i, j:end_j] = cosine_similarity(embeddings[i:end_i], embeddings[j:end_j])
    return sim


def greedy_clustering_from_matrix(sim_matrix, threshold):
    n = len(sim_matrix)
    sorted_indices = np.argsort(-np.linalg.norm(sim_matrix, axis=1))
    clusters, cluster_labels, cluster_id = {}, -1 * np.ones(n, dtype=int), 0
    for i in tqdm(sorted_indices, desc="Greedy clustering", leave=False):
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


def cluster_full_set(embeddings, threshold, cluster_path):
    """Greedy-cluster the full SAM set once; cache for reuse."""
    if os.path.exists(cluster_path):
        print(f"Loading cached clusters from {cluster_path}")
        with open(cluster_path, "rb") as f:
            return pickle.load(f)["clusters"]
    sim = batched_cosine_similarity(embeddings, save_path=str(OUT_DIR / "full_cosine_sim.dat"))
    clusters, cluster_labels = greedy_clustering_from_matrix(sim, threshold)
    with open(cluster_path, "wb") as f:
        pickle.dump({"clusters": clusters, "cluster_labels": cluster_labels}, f)
    return clusters


# ============================================================================
# Splitting
# ============================================================================

def split_clusters_train_test(clusters, y_labels, train_ratio=0.9, max_attempts=1000):
    """Split clusters into train / test, ensuring every label has a positive in test."""
    cluster_ids = list(clusters.keys())
    n_labels = y_labels.shape[1]
    for attempt in range(max_attempts):
        random.shuffle(cluster_ids)
        k = int(train_ratio * len(cluster_ids))
        train_idx = [i for c in cluster_ids[:k] for i in clusters[c]]
        test_idx = [i for c in cluster_ids[k:] for i in clusters[c]]
        if all(y_labels[test_idx][:, j].sum() > 0 for j in range(n_labels)):
            return train_idx, test_idx
    raise ValueError("Could not find a valid cluster split.")


def split_random_train_test(n_samples, y_labels, train_ratio=0.9, max_attempts=1000):
    """Random train / test split, ensuring every label has a positive in test."""
    n_labels = y_labels.shape[1]
    all_idx = list(range(n_samples))
    for attempt in range(max_attempts):
        random.shuffle(all_idx)
        k = int(train_ratio * n_samples)
        train_idx, test_idx = all_idx[:k], all_idx[k:]
        if all(y_labels[test_idx][:, j].sum() > 0 for j in range(n_labels)):
            return train_idx, test_idx
    raise ValueError("Could not find a valid random split.")


def split_train_pool_for_epoch(n_samples, train_ratio=0.9):
    idx = list(range(n_samples))
    random.shuffle(idx)
    k = int(train_ratio * n_samples)
    return idx[:k], idx[k:]


# ============================================================================
# Metrics
# ============================================================================

def compute_auprc(y_true, probs):
    vals = [auc(*precision_recall_curve(y_true[:, i], probs[:, i])[1::-1])
            for i in range(y_true.shape[1]) if len(np.unique(y_true[:, i])) > 1]
    return np.mean(vals) if vals else 0.0


def compute_auroc(y_true, probs):
    vals = [roc_auc_score(y_true[:, i], probs[:, i])
            for i in range(y_true.shape[1]) if len(np.unique(y_true[:, i])) > 1]
    return np.mean(vals) if vals else 0.0


def compute_mcc(y_true, y_pred):
    vals = [matthews_corrcoef(y_true[:, i], y_pred[:, i])
            for i in range(y_true.shape[1]) if len(np.unique(y_true[:, i])) > 1]
    return np.mean(vals) if vals else 0.0


def compute_f1(y_true, y_pred):
    vals = [f1_score(y_true[:, i], y_pred[:, i])
            for i in range(y_true.shape[1]) if len(np.unique(y_true[:, i])) > 1]
    return np.mean(vals) if vals else 0.0


def compute_accuracy(y_true, y_pred):
    return np.mean([accuracy_score(y_true[:, i], y_pred[:, i]) for i in range(y_true.shape[1])])


def compute_all_metrics(y_true, probs, threshold=0.5):
    y_pred = (probs >= threshold).astype(int)
    return {"auprc": compute_auprc(y_true, probs), "mcc": compute_mcc(y_true, y_pred),
            "auroc": compute_auroc(y_true, probs), "f1": compute_f1(y_true, y_pred),
            "accuracy": compute_accuracy(y_true, y_pred)}


def compute_per_label_metrics(y_true, probs, labels, threshold=0.5):
    y_pred = (probs >= threshold).astype(int)
    out = {}
    for i, label in enumerate(labels):
        if len(np.unique(y_true[:, i])) > 1:
            p, r, _ = precision_recall_curve(y_true[:, i], probs[:, i])
            out[label] = {"auprc": auc(r, p), "auroc": roc_auc_score(y_true[:, i], probs[:, i]),
                          "mcc": matthews_corrcoef(y_true[:, i], y_pred[:, i]),
                          "f1": f1_score(y_true[:, i], y_pred[:, i]),
                          "accuracy": accuracy_score(y_true[:, i], y_pred[:, i])}
        else:
            out[label] = {k: np.nan for k in ["auprc", "auroc", "mcc", "f1", "accuracy"]}
    return out


def find_optimal_thresholds(y_true, probs, labels):
    """Per-class MCC-optimal threshold, tuned on TRAIN predictions.

    Matches the per-label MCC-optimal thresholding used by the deployment
    pipelines and the manuscript's headline metric (MCC). Tuning on TRAIN
    (not the held-out test) avoids eval-set leakage."""
    thresholds = np.full(len(labels), 0.5)
    for i in range(len(labels)):
        scores, true = probs[:, i], y_true[:, i]
        if len(np.unique(true)) < 2:
            continue
        best_mcc, best_t = -2.0, 0.5
        for t in np.unique(np.percentile(scores, np.linspace(0, 100, 200))):
            m = matthews_corrcoef(true, (scores >= t).astype(int))
            if m > best_mcc:
                best_mcc, best_t = m, t
        thresholds[i] = best_t
    return thresholds


# ============================================================================
# MLP
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
        return torch.sigmoid(x) if self.use_sigmoid else x


def compute_pos_weight(y_train):
    counts = y_train.sum(axis=0).clip(min=1)
    return torch.tensor((len(y_train) - counts) / counts, dtype=torch.float32)


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
                weight = ti * self.pos_weight[i] + (1 - ti)
                total += (((1 - p_t) ** self.gamma) * bce * weight).mean()
            elif lt == "sqrt_weighted_bce":
                total += nn.functional.binary_cross_entropy_with_logits(
                    li, ti, pos_weight=torch.sqrt(self.pos_weight[i:i + 1]))
            elif lt == "bce":
                total += nn.functional.binary_cross_entropy_with_logits(li, ti)
            else:  # weighted_bce
                total += nn.functional.binary_cross_entropy_with_logits(
                    li, ti, pos_weight=self.pos_weight[i:i + 1])
        return total


def train_mlp(X_train, y_train, label_loss_map, epochs, lr):
    """Train the MLP on-device with per-epoch 90% resampling; return the model.

    Tensors are moved to DEVICE once and resampled with an on-device permutation
    (no per-epoch host<->device copies), so this runs on GPU when available."""
    torch.set_num_threads(20)
    model = MultilabelNN(X_train.shape[1], y_train.shape[1], use_sigmoid=False).to(DEVICE)
    criterion = PerLabelMixedLoss(
        compute_pos_weight(y_train).to(DEVICE), label_loss_map, gamma=2.0)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    Xtr = torch.tensor(X_train, dtype=torch.float32, device=DEVICE)
    ytr = torch.tensor(y_train, dtype=torch.float32, device=DEVICE)
    n = X_train.shape[0]
    k = max(1, int(0.9 * n))
    for _ in tqdm(range(epochs), desc="  MLP", leave=False, mininterval=2.0):
        idx = torch.randperm(n, device=DEVICE)[:k]
        model.train()
        optimizer.zero_grad()
        criterion(model(Xtr[idx]), ytr[idx]).backward()
        optimizer.step()
    return model


def mlp_predict(model, X):
    model.eval()
    with torch.no_grad():
        xb = torch.tensor(X, dtype=torch.float32, device=DEVICE)
        return torch.sigmoid(model(xb)).cpu().numpy()


# ============================================================================
# sklearn classifiers
# ============================================================================

def train_sklearn(clf_name, clf_template, X_train, y_train, X_test):
    """Fit a sklearn classifier; return (train_probs, test_probs)."""
    n_labels = y_train.shape[1]
    if clf_name.startswith("SVM"):
        def train_single(j):
            if len(np.unique(y_train[:, j])) < 2:
                return None
            m = clone(clf_template)
            m.fit(X_train, y_train[:, j])
            return m
        trained = Parallel(n_jobs=-1)(delayed(train_single)(j) for j in range(n_labels))
        train_probs = np.zeros((X_train.shape[0], n_labels))
        test_probs = np.zeros((X_test.shape[0], n_labels))
        for j, m in enumerate(trained):
            if m is not None:
                train_probs[:, j] = m.predict_proba(X_train)[:, 1]
                test_probs[:, j] = m.predict_proba(X_test)[:, 1]
    else:
        model = MultiOutputClassifier(clone(clf_template))
        model.fit(X_train, y_train)
        train_probs = np.hstack([e.predict_proba(X_train)[:, 1].reshape(-1, 1)
                                 for e in model.estimators_])
        test_probs = np.hstack([e.predict_proba(X_test)[:, 1].reshape(-1, 1)
                                for e in model.estimators_])
    return train_probs, test_probs


# ============================================================================
# MMseqs2 1-NN
# ============================================================================

def _mmseqs_vote(result_tsv, query_to_idx, n_query, train_label_map, n_labels, exclude_self=False):
    y = np.zeros((n_query, n_labels), dtype=np.float32)
    if os.path.getsize(result_tsv) == 0:
        return y
    res = pd.read_csv(result_tsv, sep="\t", header=None, names=["query", "target", "evalue", "bits"])
    for query, group in res.groupby("query"):
        query = str(query)
        if query not in query_to_idx:
            continue
        idx = query_to_idx[query]
        if exclude_self:
            group = group[group["target"].astype(str) != query]
        if group.empty:
            continue
        top = group.nlargest(10, "bits")
        total = top["bits"].values.clip(min=0).sum()
        if total > 0:
            for _, row in top.iterrows():
                tgt = str(row["target"])
                if tgt in train_label_map and row["bits"] > 0:
                    y[idx] += (row["bits"] / total) * train_label_map[tgt]
        else:
            tgt = str(top.iloc[0]["target"])
            if tgt in train_label_map:
                y[idx] = train_label_map[tgt]
    return y


def run_mmseqs2(train_entries, train_seqs, train_labels, test_entries, test_seqs, test_labels, labels):
    n_test = len(test_entries)
    with tempfile.TemporaryDirectory() as tmp:
        train_fa, test_fa = os.path.join(tmp, "train.fasta"), os.path.join(tmp, "test.fasta")
        with open(train_fa, "w") as f:
            for e, s in zip(train_entries, train_seqs):
                f.write(f">{e}\n{s}\n")
        with open(test_fa, "w") as f:
            for e, s in zip(test_entries, test_seqs):
                f.write(f">{e}\n{s}\n")
        trainDB, testDB = os.path.join(tmp, "trainDB"), os.path.join(tmp, "testDB")
        resDB, tdir = os.path.join(tmp, "resDB"), os.path.join(tmp, "tmp")
        subprocess.run(["mmseqs", "createdb", train_fa, trainDB], capture_output=True, check=True)
        subprocess.run(["mmseqs", "createdb", test_fa, testDB], capture_output=True, check=True)
        subprocess.run(["mmseqs", "search", testDB, trainDB, resDB, tdir,
                        "--max-seqs", "10", "-e", "100", "--threads", "20"],
                       capture_output=True, check=True)
        res_tsv = os.path.join(tmp, "res.tsv")
        subprocess.run(["mmseqs", "convertalis", testDB, trainDB, resDB, res_tsv,
                        "--format-output", "query,target,evalue,bits"], capture_output=True, check=True)
        # train-vs-train (self excluded) for threshold tuning
        resTrainDB, res_train_tsv = os.path.join(tmp, "resTrainDB"), os.path.join(tmp, "res_train.tsv")
        subprocess.run(["mmseqs", "search", trainDB, trainDB, resTrainDB, os.path.join(tmp, "tmp2"),
                        "--max-seqs", "11", "-e", "100", "--threads", "20"], capture_output=True, check=True)
        subprocess.run(["mmseqs", "convertalis", trainDB, trainDB, resTrainDB, res_train_tsv,
                        "--format-output", "query,target,evalue,bits"], capture_output=True, check=True)

        train_label_map = {e: row for e, row in zip(train_entries, train_labels)}
        test_to_idx = {e: i for i, e in enumerate(test_entries)}
        train_to_idx = {e: i for i, e in enumerate(train_entries)}
        test_probs = _mmseqs_vote(res_tsv, test_to_idx, n_test, train_label_map, len(labels))
        train_probs = _mmseqs_vote(res_train_tsv, train_to_idx, len(train_entries),
                                   train_label_map, len(labels), exclude_self=True)

    thresholds = find_optimal_thresholds(train_labels, train_probs, labels)
    return (test_probs, compute_all_metrics(test_labels, test_probs, thresholds),
            compute_per_label_metrics(test_labels, test_probs, labels, thresholds))


# ============================================================================
# Foldseek 1-NN (structural-homology baseline; structural twin of MMseqs2 1-NN)
# ============================================================================

def _name_to_entry(name):
    """Foldseek query/target name (e.g. 'P12345.pdb_A') -> UniProt Entry.

    Filenames are {Entry}.pdb and UniProt accessions contain no '_', so
    stripping '.pdb' and any trailing '_<chain>' recovers the Entry."""
    name = str(name).split()[0].replace(".pdb", "")
    return name.split("_")[0]


def _foldseek_db(entries, struct_dir, dbpath, tmp):
    """Symlink the needed {Entry}.pdb into a dir and foldseek-createdb it."""
    sdir = os.path.join(tmp, os.path.basename(dbpath) + "_pdbs")
    os.makedirs(sdir, exist_ok=True)
    n = 0
    for e in entries:
        src = os.path.join(struct_dir, f"{e}.pdb")
        if os.path.exists(src):
            os.symlink(os.path.abspath(src), os.path.join(sdir, f"{e}.pdb"))
            n += 1
    subprocess.run([FOLDSEEK, "createdb", sdir, dbpath], capture_output=True, check=True)
    return n


def _foldseek_vote(res, query_to_idx, n_query, train_label_map, n_labels, exclude_self=False):
    """Top-10 bit-score-weighted label vote (identical scheme to _mmseqs_vote),
    over a results frame whose query/target are already mapped to Entries."""
    y = np.zeros((n_query, n_labels), dtype=np.float32)
    if res is None or res.empty:
        return y
    for query, group in res.groupby("query"):
        if query not in query_to_idx:
            continue
        idx = query_to_idx[query]
        if exclude_self:
            group = group[group["target"] != query]
        if group.empty:
            continue
        top = group.nlargest(10, "bits")
        total = top["bits"].values.clip(min=0).sum()
        if total > 0:
            for _, row in top.iterrows():
                if row["target"] in train_label_map and row["bits"] > 0:
                    y[idx] += (row["bits"] / total) * train_label_map[row["target"]]
        else:
            t0 = top.iloc[0]["target"]
            if t0 in train_label_map:
                y[idx] = train_label_map[t0]
    return y


def run_foldseek(train_entries, train_labels, test_entries, test_labels, labels, struct_dir):
    """Structural 1-NN: vote substrate type from a query's nearest structural
    neighbours in the train set (Foldseek 3Di alignment). Thresholds tuned on
    train-vs-train (self excluded), mirroring run_mmseqs2."""
    n_test = len(test_entries)
    with tempfile.TemporaryDirectory() as tmp:
        trainDB, testDB = os.path.join(tmp, "trainDB"), os.path.join(tmp, "testDB")
        _foldseek_db(train_entries, struct_dir, trainDB, tmp)
        _foldseek_db(test_entries, struct_dir, testDB, tmp)

        resDB, res_tsv = os.path.join(tmp, "resDB"), os.path.join(tmp, "res.tsv")
        subprocess.run([FOLDSEEK, "search", testDB, trainDB, resDB, os.path.join(tmp, "t1"),
                        "-e", "100", "--max-seqs", "10", "--threads", "16"],
                       capture_output=True, check=True)
        subprocess.run([FOLDSEEK, "convertalis", testDB, trainDB, resDB, res_tsv,
                        "--format-output", "query,target,evalue,bits"],
                       capture_output=True, check=True)
        # train-vs-train (self excluded) for threshold tuning
        resTrainDB, res_train_tsv = os.path.join(tmp, "resTrainDB"), os.path.join(tmp, "res_train.tsv")
        subprocess.run([FOLDSEEK, "search", trainDB, trainDB, resTrainDB, os.path.join(tmp, "t2"),
                        "-e", "100", "--max-seqs", "11", "--threads", "16"],
                       capture_output=True, check=True)
        subprocess.run([FOLDSEEK, "convertalis", trainDB, trainDB, resTrainDB, res_train_tsv,
                        "--format-output", "query,target,evalue,bits"],
                       capture_output=True, check=True)

        def load(tsv):
            if os.path.getsize(tsv) == 0:
                return None
            df = pd.read_csv(tsv, sep="\t", header=None, names=["query", "target", "evalue", "bits"])
            df["query"] = df["query"].map(_name_to_entry)
            df["target"] = df["target"].map(_name_to_entry)
            return df

        res, res_train = load(res_tsv), load(res_train_tsv)
        train_label_map = {e: row for e, row in zip(train_entries, train_labels)}
        test_to_idx = {e: i for i, e in enumerate(test_entries)}
        train_to_idx = {e: i for i, e in enumerate(train_entries)}
        test_probs = _foldseek_vote(res, test_to_idx, n_test, train_label_map, len(labels))
        train_probs = _foldseek_vote(res_train, train_to_idx, len(train_entries),
                                     train_label_map, len(labels), exclude_self=True)

    thresholds = find_optimal_thresholds(train_labels, train_probs, labels)
    return (test_probs, compute_all_metrics(test_labels, test_probs, thresholds),
            compute_per_label_metrics(test_labels, test_probs, labels, thresholds))


# ============================================================================
# Profile HMM
# ============================================================================

def _parse_tblout(path, entry_to_idx, scores, class_idx):
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            target, score = parts[0], float(parts[5])
            if target in entry_to_idx:
                idx = entry_to_idx[target]
                scores[idx, class_idx] = max(scores[idx, class_idx], score)


def run_profile_hmm(train_entries, train_seqs, train_labels, test_entries, test_seqs,
                    test_labels, labels, max_hmm_seqs=500):
    n_test, n_train = len(test_entries), len(train_entries)
    with tempfile.TemporaryDirectory() as tmp:
        test_fa, train_fa = os.path.join(tmp, "test.fasta"), os.path.join(tmp, "train.fasta")
        with open(test_fa, "w") as f:
            for e, s in zip(test_entries, test_seqs):
                f.write(f">{e}\n{s}\n")
        with open(train_fa, "w") as f:
            for e, s in zip(train_entries, train_seqs):
                f.write(f">{e}\n{s}\n")
        test_to_idx = {e: i for i, e in enumerate(test_entries)}
        train_to_idx = {e: i for i, e in enumerate(train_entries)}
        test_scores = np.full((n_test, len(labels)), -1000.0)
        train_scores = np.full((n_train, len(labels)), -1000.0)

        for class_idx, label in enumerate(labels):
            mask = train_labels[:, class_idx] == 1
            pos_e = [e for e, m in zip(train_entries, mask) if m]
            pos_s = [s for s, m in zip(train_seqs, mask) if m]
            if len(pos_e) < 2:
                continue
            if len(pos_e) > max_hmm_seqs:
                samp = random.sample(range(len(pos_e)), max_hmm_seqs)
                pos_e, pos_s = [pos_e[i] for i in samp], [pos_s[i] for i in samp]
            class_fa = os.path.join(tmp, f"{label}.fasta")
            with open(class_fa, "w") as f:
                for e, s in zip(pos_e, pos_s):
                    f.write(f">{e}\n{s}\n")
            aligned = os.path.join(tmp, f"{label}_aln.fasta")
            if subprocess.run(f"mafft --auto --thread 20 --quiet {class_fa} > {aligned}",
                              shell=True, capture_output=True).returncode != 0:
                continue
            hmm = os.path.join(tmp, f"{label}.hmm")
            if subprocess.run(["hmmbuild", "--amino", "--cpu", "20", hmm, aligned],
                              capture_output=True).returncode != 0:
                continue
            test_tbl = os.path.join(tmp, f"{label}_test.tbl")
            subprocess.run(["hmmsearch", "--tblout", test_tbl, "-E", "1000", "--cpu", "20", hmm, test_fa],
                           capture_output=True, check=True)
            _parse_tblout(test_tbl, test_to_idx, test_scores, class_idx)
            train_tbl = os.path.join(tmp, f"{label}_train.tbl")
            subprocess.run(["hmmsearch", "--tblout", train_tbl, "-E", "1000", "--cpu", "20", hmm, train_fa],
                           capture_output=True, check=True)
            _parse_tblout(train_tbl, train_to_idx, train_scores, class_idx)

        thresholds = np.zeros(len(labels))
        for class_idx in range(len(labels)):
            scores, true = train_scores[:, class_idx], train_labels[:, class_idx]
            valid = scores > -999
            if valid.sum() == 0:
                continue
            best_mcc, best_t = -2.0, 0
            for t in np.percentile(scores[valid], np.linspace(0, 100, 200)):
                m = matthews_corrcoef(true, (scores >= t).astype(int))
                if m > best_mcc:
                    best_mcc, best_t = m, t
            thresholds[class_idx] = best_t

        test_probs = np.zeros((n_test, len(labels)), dtype=np.float32)
        for ci in range(len(labels)):
            valid_train = train_scores[:, ci][train_scores[:, ci] > -999]
            scale = np.std(valid_train) * 2 if len(valid_train) > 1 and np.std(valid_train) > 0 else 50.0
            test_probs[:, ci] = 1 / (1 + np.exp(-(test_scores[:, ci] - thresholds[ci]) / scale))
        y_pred = np.zeros((n_test, len(labels)), dtype=int)
        for ci in range(len(labels)):
            y_pred[:, ci] = (test_scores[:, ci] >= thresholds[ci]).astype(int)

    metrics = {"auprc": compute_auprc(test_labels, test_probs),
               "auroc": compute_auroc(test_labels, test_probs),
               "mcc": compute_mcc(test_labels, y_pred), "f1": compute_f1(test_labels, y_pred),
               "accuracy": compute_accuracy(test_labels, y_pred)}
    per_label = {}
    for i, label in enumerate(labels):
        if len(np.unique(test_labels[:, i])) > 1:
            p, r, _ = precision_recall_curve(test_labels[:, i], test_probs[:, i])
            per_label[label] = {"auprc": auc(r, p), "auroc": roc_auc_score(test_labels[:, i], test_probs[:, i]),
                                "mcc": matthews_corrcoef(test_labels[:, i], y_pred[:, i]),
                                "f1": f1_score(test_labels[:, i], y_pred[:, i]),
                                "accuracy": accuracy_score(test_labels[:, i], y_pred[:, i])}
        else:
            per_label[label] = {k: np.nan for k in ["auprc", "auroc", "mcc", "f1", "accuracy"]}
    return test_probs, metrics, per_label


# ============================================================================
# One dataset x split x iteration: train all classifiers, score on shared test
# ============================================================================

def replay_hmm_rng(y_train, max_hmm_seqs=500):
    """Consume exactly the random draws run_profile_hmm() makes for this training set (one
    random.sample per label with more than max_hmm_seqs positives), so a run that skips the
    HMM still produces the same splits for all later iterations."""
    for class_idx in range(y_train.shape[1]):
        n_pos = int((y_train[:, class_idx] == 1).sum())
        if n_pos >= 2 and n_pos > max_hmm_seqs:
            random.sample(range(n_pos), max_hmm_seqs)


def evaluate_all_classifiers(dataset, X_train, y_train, entries_train, seqs_train,
                             X_test, y_test, entries_test, seqs_test,
                             src_train, src_test):
    """Train every classifier on the given train set and score on the shared test set.
    src_train / src_test are {embedding: matrix} dicts (row-aligned to train / test)
    used by the MLP family. Returns (overall, per_label) dicts keyed by classifier."""
    overall, per_label = {}, {}

    # sklearn (probabilistic): tune thresholds on train predictions
    for clf_name, template in ({} if ONLY_MLP else SKLEARN_CLASSIFIERS).items():
        print(f"    [{dataset}] {clf_name}")
        train_probs, test_probs = train_sklearn(clf_name, template, X_train, y_train, X_test)
        thr = find_optimal_thresholds(y_train, train_probs, LABELS)
        overall[clf_name] = compute_all_metrics(y_test, test_probs, thr)
        per_label[clf_name] = compute_per_label_metrics(y_test, test_probs, LABELS, thr)

    # MLP family: one per embedding, each with its own per-dataset epoch count
    # and per-label mixed loss. Thresholds tuned on train predictions.
    for mlp_name, cfg in MLP_VARIANTS.items():
        dcfg = cfg[dataset]
        Xtr = np.hstack([src_train[s] for s in cfg["sources"]])
        Xte = np.hstack([src_test[s] for s in cfg["sources"]])
        print(f"    [{dataset}] {mlp_name}  (dim={Xtr.shape[1]}, {dcfg['epochs']} ep)")
        model = train_mlp(Xtr, y_train, dcfg["loss"], dcfg["epochs"], LEARNING_RATE)
        train_probs, test_probs = mlp_predict(model, Xtr), mlp_predict(model, Xte)
        thr = find_optimal_thresholds(y_train, train_probs, LABELS)
        overall[mlp_name] = compute_all_metrics(y_test, test_probs, thr)
        per_label[mlp_name] = compute_per_label_metrics(y_test, test_probs, LABELS, thr)

    if ONLY_MLP:
        replay_hmm_rng(y_train)          # keep the random stream (and so the splits) identical
        return overall, per_label

    # MMseqs2 1-NN (tunes thresholds internally on train)
    print(f"    [{dataset}] MMseqs2 1-NN")
    _, m, pl = run_mmseqs2(entries_train, seqs_train, y_train, entries_test, seqs_test, y_test, LABELS)
    overall["MMseqs2_1NN"], per_label["MMseqs2_1NN"] = m, pl

    # Profile HMM (tunes thresholds internally on train)
    print(f"    [{dataset}] Profile HMM")
    _, m, pl = run_profile_hmm(entries_train, seqs_train, y_train, entries_test, seqs_test, y_test, LABELS)
    overall["Profile_HMM"], per_label["Profile_HMM"] = m, pl

    # Foldseek 1-NN (structural-homology baseline; tunes thresholds on train)
    print(f"    [{dataset}] Foldseek 1-NN")
    _, m, pl = run_foldseek(entries_train, y_train, entries_test, y_test, LABELS, str(STRUCT_DIR))
    overall["Foldseek_1NN"], per_label["Foldseek_1NN"] = m, pl

    return overall, per_label


# ============================================================================
# Main
# ============================================================================

def main():
    if not ONLY_MLP:
        check_tools()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print("SAM classifier comparison")
    print("=" * 70)
    print(f"output: {OUT_DIR} | datasets: {DATASETS} | only MLPs: {ONLY_MLP}")
    for _n, _c in MLP_VARIANTS.items():
        _txt = " | ".join(f"{d} {_c[d]['epochs']:>4} ep {list(_c[d]['loss'].values())}" for d in DATASETS)
        print(f"  {_n:20s} {_txt}")

    print(f"Device: {DEVICE}. Loading embedding sources: {', '.join(EMB_SOURCES)}")
    source_maps = {s: load_source_map(p, layout)
                   for s, (p, layout, _dim) in EMB_SOURCES.items()}
    full_df, full_mats = attach_sources(pd.read_csv(FULL_CSV, low_memory=False), source_maps)
    bal_df, bal_mats = attach_sources(pd.read_csv(BAL_CSV, low_memory=False), source_maps)
    print(f"Full SAM: {len(full_df)} | Balanced SAM: {len(bal_df)} "
          f"(entries covered by all {len(EMB_SOURCES)} sources)")

    full_X = full_mats[SKLEARN_SOURCE]        # sklearn features
    full_y = full_df[LABELS].values.astype(np.float32)
    full_entries = full_df["Entry"].values
    full_seqs = full_df["Sequence"].values

    bal_X = bal_mats[SKLEARN_SOURCE]
    bal_y = bal_df[LABELS].values.astype(np.float32)
    bal_entries = bal_df["Entry"].values
    bal_seqs = bal_df["Sequence"].values

    print(f"sklearn features: {SKLEARN_SOURCE} (dim={full_X.shape[1]})  |  "
          f"clustering space: {CLUSTER_SOURCE} @ cosine {THRESHOLD}")
    clusters = cluster_full_set(full_mats[CLUSTER_SOURCE], THRESHOLD,
                                str(OUT_DIR / f"full_clusters_t{THRESHOLD}.pkl"))
    print(f"Full-set clusters: {len(clusters)}")

    detailed, per_label_rows = [], []

    for split_mode in SPLIT_MODES:
        for iteration in range(N_ITERATIONS):
            print("\n" + "=" * 70)
            print(f"SPLIT={split_mode.upper()}  ITERATION {iteration + 1}/{N_ITERATIONS}")
            print("=" * 70)

            # 1. test set from the full-SAM split
            if split_mode == "cluster":
                train_idx, test_idx = split_clusters_train_test(clusters, full_y, TRAIN_RATIO)
            else:
                train_idx, test_idx = split_random_train_test(len(full_X), full_y, TRAIN_RATIO)
            test_entry_set = set(full_entries[test_idx])
            print(f"  full train={len(train_idx)}  shared test={len(test_idx)}")

            X_test = full_X[test_idx]
            y_test = full_y[test_idx]
            e_test = full_entries[test_idx]
            s_test = full_seqs[test_idx]
            src_test = {s: full_mats[s][test_idx] for s in EMB_SOURCES}

            # 2. full training set
            full_train = {
                "X": full_X[train_idx], "y": full_y[train_idx],
                "e": full_entries[train_idx], "s": full_seqs[train_idx],
                "src": {s: full_mats[s][train_idx] for s in EMB_SOURCES},
            }

            # 3. balanced training set = balanced minus any test-set entries
            keep = np.array([e not in test_entry_set for e in bal_entries])
            bal_train = {"X": bal_X[keep], "y": bal_y[keep],
                         "e": bal_entries[keep], "s": bal_seqs[keep],
                         "src": {s: bal_mats[s][keep] for s in EMB_SOURCES}}
            print(f"  balanced train={keep.sum()} (withheld {len(bal_entries) - keep.sum()} test entries)")

            for dataset, tr in [("full", full_train), ("balanced", bal_train)]:
                if dataset not in DATASETS:
                    replay_hmm_rng(tr["y"])      # a skipped arm still advances the random stream, so later splits are unchanged
                    continue
                overall, per_label = evaluate_all_classifiers(
                    dataset, tr["X"], tr["y"], tr["e"], tr["s"],
                    X_test, y_test, e_test, s_test, tr["src"], src_test)
                for clf_name, metrics in overall.items():
                    detailed.append({"split_mode": split_mode, "dataset": dataset,
                                     "iteration": iteration + 1, "classifier": clf_name, **metrics})
                    for label, lm in per_label[clf_name].items():
                        per_label_rows.append({"split_mode": split_mode, "dataset": dataset,
                                               "iteration": iteration + 1, "classifier": clf_name,
                                               "label": label, **lm})

    # ---- Save ----
    detailed_df = pd.DataFrame(detailed)
    detailed_df.to_csv(OUT_DIR / "sam_classifier_comparison_detailed.csv", index=False)
    pd.DataFrame(per_label_rows).to_csv(OUT_DIR / "sam_classifier_comparison_per_label.csv", index=False)

    metric_cols = ["auprc", "auroc", "mcc", "f1", "accuracy"]
    summary = (detailed_df.groupby(["split_mode", "dataset", "classifier"])[metric_cols]
               .agg(["mean", "std"]).reset_index())
    summary.columns = ["_".join(c).strip("_") for c in summary.columns]
    summary.to_csv(OUT_DIR / "sam_classifier_comparison_summary.csv", index=False)
    print(f"\nWrote detailed / summary / per-label CSVs to {OUT_DIR}")

    # ---- Figure: mean metric bars (rows=metric, cols=split x dataset) ----
    panels = [(sm, ds) for sm in SPLIT_MODES for ds in DATASETS]
    fig, axes = plt.subplots(len(metric_cols), len(panels),
                             figsize=(5 * len(panels), 3.2 * len(metric_cols)), squeeze=False)
    for r, metric in enumerate(metric_cols):
        for c, (sm, ds) in enumerate(panels):
            ax = axes[r][c]
            sub = detailed_df[(detailed_df.split_mode == sm) & (detailed_df.dataset == ds)]
            g = sub.groupby("classifier")[metric].agg(["mean", "std"]).reindex(
                list(SKLEARN_CLASSIFIERS) + list(MLP_VARIANTS)
                + ["MMseqs2_1NN", "Foldseek_1NN", "Profile_HMM"])
            ax.bar(range(len(g)), g["mean"], yerr=g["std"], capsize=3)
            ax.set_xticks(range(len(g)))
            ax.set_xticklabels(g.index, rotation=90, fontsize=7)
            ax.set_ylim(0, 1)
            ax.set_ylabel(metric.upper())
            if r == 0:
                ax.set_title(f"{sm} / {ds}")
            ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "sam_classifier_comparison.png", dpi=150)
    plt.close()
    print(f"Wrote {OUT_DIR / 'sam_classifier_comparison.png'}")
    print("\n" + "=" * 70 + "\nCOMPLETE\n" + "=" * 70)


if __name__ == "__main__":
    main()
