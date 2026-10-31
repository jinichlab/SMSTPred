#!/usr/bin/env python3
"""
Train vs. VALIDATION loss per epoch -- the overfitting counterpart to the
MCC-vs-epoch grids of Figures S2 (SAM full) and S5 (UP ProstT5); the output is Figures S3 and S6.
(The balanced SAM subset can be run too, tasks 6-11; it is not used in the manuscript.)

The MCC grids are scored only at checkpoints and only as a metric. This script
re-runs the *same* sweeps (same embeddings, same loss types, same 6
cluster-disjoint splits, same MLP, same 20,000-epoch budget) but also records,
after every logged epoch, the loss on the full training pool and on the
cluster-disjoint validation set, per label. Overfitting shows up as validation
loss turning upward while training loss keeps falling.

Why a rerun and not the old in-loop "test loss": the earlier scripts'
`test_losses` are a random 10% slice of the TRAIN pool, redrawn every epoch (so
the same protein is trained on one step and "tested" the next) and not
cluster-disjoint. They cannot diverge honestly. Here the validation loss is on
the held-out clusters only.

Protocol reproduced exactly from the sweep scripts
  SAM (S2):    mlp_embedding_comparison.py (5 embeddings) and mlp_embedding_comparison_saprot.py (SaProt). Cluster ONCE on
               ProtT5 @ cosine 0.6 over each dataset's common-intersection set,
               6 cluster-disjoint 90/10 splits drawn from the seed-42 stream
               (full first, then balanced -- the balanced splits depend on the
               RNG state left by the full ones, so full is always regenerated).
  UP  (S5):    mlp_loss_comparison_prostt5.py. ProstT5 on UP, cached ProtT5 @ 0.5
               cluster ids, 6 splits from the seed-42 stream.
  Model/loss:  MLP in -> 512 -> 256 -> 4 (raw logits), Adam lr 1e-3, ONE gradient
               step per epoch on a fresh random 90% of the train pool, loss =
               PerLabelUniformLoss (bce / weighted_bce / sqrt_weighted_bce /
               focal), summed over the 4 labels. Here it is returned per label so
               the label losses can be shown separately; their sum is the optimized
               objective (checked equal to the original class).

Losses are evaluated with the same criterion and the same TRAIN-pool pos_weight
on both sets, so train and validation numbers share a scale. Validation loss is
NOT used to pick anything here.

Sanity check built in: at the original checkpoint epochs the script also records
val per-label MCC (SAM: MCC-optimal thresholds tuned on TRAIN, as in the sweep;
UP: fixed 0.5), so the rerun can be compared against the S2/S5 CSVs. Seeds
differ from the original run (the sweep's RNG stream is not reproducible once
tasks are parallelised), so agreement is within seed noise, not exact.

Tasks (SLURM array, see run_train_val_loss.sh):  0-5 SAM full x 6 embeddings,
6-11 SAM balanced x 6 embeddings, 12-15 UP ProstT5 x 4 losses.

Env: smstpred (`make env`) (GPU).
Outputs (./train_val_loss_output/):
  <set>__<embedding>[__<loss>]__curves.csv.gz   dataset, embedding, loss_type, iteration, epoch, label, train_loss, val_loss
  <set>__<embedding>[__<loss>]__ckpt_mcc.csv.gz dataset, embedding, loss_type, iteration, epoch, label, mcc
"""
import argparse
import os
import random
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics.pairwise import cosine_similarity

# ============================================================================
# Configuration
# ============================================================================
# Inputs: the Zenodo files, in <repo>/data or $SMSTPRED_DATA (see data/README.md).
FILES_DIR = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[1] / "data"))
UP_CLUSTER_CSV = FILES_DIR / "clusters_UP_ProtT5_t0.5.csv"       # Entry, cluster (row-aligned to labels_UP.csv)
OUT_DIR = Path(__file__).resolve().parent / "train_val_loss_output"

SAM_CSV = {
    "full":     FILES_DIR / "labels_SAM_full.csv",        # deduped: 7,405 proteins
    "balanced": FILES_DIR / "labels_SAM_balanced.csv",    # deduped: 2,924 proteins
}
UP_CSV = FILES_DIR / "labels_UP.csv"
UP_PROSTT5_H5 = FILES_DIR / "embeddings_UP_ProstT5.h5"

# group -> embedding sources ((path, h5 layout)) and the conditions built from them.
# Two groups mirror the two original scripts: their common-intersection sets differ
# (ESM-IF drops ~26 proteins without an AFDB structure; SaProt covers all 7,405), so
# their splits are NOT interchangeable.
SAM_GROUPS = {
    "main": {
        "sources": {"prott5":  (FILES_DIR / "embeddings_SAM_ProtT5.h5",         "flat"),
                    "prostt5": (FILES_DIR / "embeddings_SAM_ProstT5.h5", "mean_group"),
                    "esmif":   (FILES_DIR / "embeddings_SAM_ESM-IF.h5",   "flat")},
        "conditions": {"ProtT5": ["prott5"], "ProstT5": ["prostt5"], "ESM-IF": ["esmif"],
                       "ProtT5+ESM-IF": ["prott5", "esmif"],
                       "ProstT5+ESM-IF": ["prostt5", "esmif"]},
    },
    "saprot": {
        "sources": {"prott5": (FILES_DIR / "embeddings_SAM_ProtT5.h5",        "flat"),   # clustering only
                    "saprot": (FILES_DIR / "embeddings_SAM_SaProt.h5", "flat")},
        "conditions": {"SaProt": ["saprot"]},
    },
}
EMBEDDINGS = ["ProstT5", "SaProt", "ProtT5", "ESM-IF", "ProtT5+ESM-IF", "ProstT5+ESM-IF"]
GROUP_OF = {c: g for g, cfg in SAM_GROUPS.items() for c in cfg["conditions"]}

LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
LOSS_TYPES = ["bce", "weighted_bce", "sqrt_weighted_bce", "focal"]
SPLIT_SEED = 42
SEED = 42                   # base for per-(iteration, loss) torch seeds
SAM_THRESHOLD = 0.6
N_ITERATIONS = 6
LEARNING_RATE = 0.001
MAX_EPOCHS = 20000
# Dense logging: every epoch to 1k (resolves the fast early dynamics), every 10 after.
LOG_EPOCHS = list(range(1, 1001)) + list(range(1010, MAX_EPOCHS + 1, 10))
# Original checkpoint grid of the S2/S5 sweeps (all members of LOG_EPOCHS).
EPOCH_CHECKPOINTS = list(range(100, 1000, 100)) + list(range(1000, 21000, 1000))

TASKS = ([("sam_full", e, None) for e in EMBEDDINGS] +
         [("sam_balanced", e, None) for e in EMBEDDINGS] +
         [("up", "ProstT5", l) for l in LOSS_TYPES])

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# ============================================================================
# Data + splits (copied from the sweep scripts so the splits are reproduced)
# ============================================================================

def load_source(path, layout):
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            item = f[k]
            if layout == "flat":
                out[k] = np.asarray(item[:], dtype=np.float32)
            elif layout == "mean_group":
                out[k] = np.asarray(item["mean"][:], dtype=np.float32)
            else:
                raise ValueError(layout)
    return out


def greedy_clustering_from_matrix(sim_matrix, threshold):
    n = len(sim_matrix)
    order = np.argsort(-np.linalg.norm(sim_matrix, axis=1))
    clusters, labels, cid = {}, -1 * np.ones(n, dtype=int), 0
    for i in order:
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
    for _ in range(max_attempts):
        random.shuffle(cids)
        k = int(train_ratio * len(cids))
        tr = [i for c in cids[:k] for i in clusters[c]]
        va = [i for c in cids[k:] for i in clusters[c]]
        if all(y[va][:, j].sum() > 0 for j in range(y.shape[1])):
            return tr, va
    raise ValueError("no valid cluster split")


def sam_dataset_and_splits(group, upto):
    """Regenerate (y, feats, splits) for dataset `upto` ('full' or 'balanced')
    exactly as the sweep script did: one seed-42 stream, full then balanced."""
    cfg = SAM_GROUPS[group]
    src = {name: load_source(p, lay) for name, (p, lay) in cfg["sources"].items()}
    random.seed(SPLIT_SEED)
    for ds in ["full", "balanced"]:
        df = pd.read_csv(SAM_CSV[ds], low_memory=False).drop_duplicates("Entry")
        label_map = {e: row for e, row in zip(df["Entry"], df[LABELS].values.astype(np.float32))}
        common = set(label_map)
        for name in src:
            common &= set(src[name])
        common = sorted(common)
        y = np.vstack([label_map[e] for e in common]).astype(np.float32)
        feats = {n: np.vstack([src[n][e] for e in common]).astype(np.float32) for n in src}
        sim = cosine_similarity(feats["prott5"])
        clusters = greedy_clustering_from_matrix(sim, SAM_THRESHOLD)
        del sim
        splits = [split_clusters_train_val(clusters, y, 0.9) for _ in range(N_ITERATIONS)]
        print(f"[{group}/{ds}] common={len(common)} clusters={len(clusters)} "
              f"val sizes={[len(v) for _, v in splits]}", flush=True)
        if ds == upto:
            return y, feats, splits


def up_dataset_and_splits():
    df = pd.read_csv(UP_CSV, low_memory=False).drop_duplicates("Entry").reset_index(drop=True)
    cl = pd.read_csv(UP_CLUSTER_CSV)
    assert list(cl["Entry"]) == list(df["Entry"]), "cluster table and label table are not row-aligned"
    entry_to_cluster = dict(zip(cl["Entry"].values, cl["cluster"].values))
    feats = load_source(UP_PROSTT5_H5, "mean_group")
    df = df[df["Entry"].apply(lambda e: e in feats)].reset_index(drop=True)
    X = np.vstack([feats[e] for e in df["Entry"]]).astype(np.float32)
    y = df[LABELS].values.astype(np.float32)
    clusters = {}
    for i, e in enumerate(df["Entry"]):
        clusters.setdefault(int(entry_to_cluster[e]), []).append(i)
    random.seed(SPLIT_SEED)
    splits = [split_clusters_train_val(clusters, y, 0.9) for _ in range(N_ITERATIONS)]
    print(f"[up] n={len(df)} clusters={len(clusters)} val sizes={[len(v) for _, v in splits]}",
          flush=True)
    return X, y, splits


# ============================================================================
# Model / loss / thresholds
# ============================================================================

class MultilabelNN(nn.Module):
    def __init__(self, input_size, output_size):
        super().__init__()
        self.fc1 = nn.Linear(input_size, 512)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, output_size)

    def forward(self, x):
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))   # raw logits


def compute_pos_weight(y):
    counts = y.sum(axis=0).clip(min=1)
    return torch.tensor((len(y) - counts) / counts, dtype=torch.float32)


class PerLabelLoss(nn.Module):
    """PerLabelUniformLoss from the sweep scripts, returning the per-label terms
    ([n_labels] tensor). .sum() is the original (optimized) objective."""
    def __init__(self, loss_type, pos_weight, gamma=2.0):
        super().__init__()
        self.loss_type, self.gamma = loss_type, gamma
        self.register_buffer("pos_weight", pos_weight)
        self.register_buffer("sqrt_pos_weight", torch.sqrt(pos_weight))

    def forward(self, logits, targets):
        terms = []
        for i in range(logits.shape[1]):
            li, ti = logits[:, i], targets[:, i]
            if self.loss_type == "bce":
                t = nn.functional.binary_cross_entropy_with_logits(li, ti)
            elif self.loss_type == "weighted_bce":
                t = nn.functional.binary_cross_entropy_with_logits(
                    li, ti, pos_weight=self.pos_weight[i:i + 1])
            elif self.loss_type == "sqrt_weighted_bce":
                t = nn.functional.binary_cross_entropy_with_logits(
                    li, ti, pos_weight=self.sqrt_pos_weight[i:i + 1])
            elif self.loss_type == "focal":
                bce = nn.functional.binary_cross_entropy_with_logits(li, ti, reduction="none")
                p = torch.sigmoid(li)
                p_t = ti * p + (1 - ti) * (1 - p)
                pw = ti * self.pos_weight[i] + (1 - ti)
                t = (((1 - p_t) ** self.gamma) * bce * pw).mean()
            terms.append(t)
        return torch.stack(terms)


def find_mcc_optimal_thresholds(y_true, probs):
    """Per-label MCC-optimal thresholds tuned on TRAIN predictions (as in the sweep)."""
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


def val_mcc_per_label(y_val, probs, thr):
    pred = (probs >= thr).astype(int)
    return [matthews_corrcoef(y_val[:, i], pred[:, i]) if len(np.unique(y_val[:, i])) > 1
            else np.nan for i in range(len(LABELS))]


# ============================================================================
# Training with per-epoch train/val loss
# ============================================================================

def train_with_loss_curves(X_tr, y_tr, X_va, y_va, loss_type, seed, tune_thr,
                           max_epochs, log_epochs, ckpt_epochs):
    """Returns (tr_hist, va_hist [n_log, n_labels], ckpt {epoch: per-label val MCC})."""
    torch.manual_seed(seed)
    model = MultilabelNN(X_tr.shape[1], len(LABELS)).to(DEVICE)
    crit = PerLabelLoss(loss_type, compute_pos_weight(y_tr)).to(DEVICE)
    opt = optim.Adam(model.parameters(), lr=LEARNING_RATE)

    Xtr = torch.tensor(X_tr, dtype=torch.float32, device=DEVICE)
    ytr = torch.tensor(y_tr, dtype=torch.float32, device=DEVICE)
    Xva = torch.tensor(X_va, dtype=torch.float32, device=DEVICE)
    yva = torch.tensor(y_va, dtype=torch.float32, device=DEVICE)
    n = Xtr.shape[0]
    k = max(1, int(0.9 * n))

    log_set, ckpt_set = set(log_epochs), set(ckpt_epochs)
    tr_hist = torch.zeros(len(log_epochs), len(LABELS), device=DEVICE)
    va_hist = torch.zeros_like(tr_hist)
    ckpt, j = {}, 0
    for epoch in range(1, max_epochs + 1):
        idx = torch.randperm(n, device=DEVICE)[:k]
        model.train(); opt.zero_grad()
        crit(model(Xtr[idx]), ytr[idx]).sum().backward(); opt.step()
        if epoch in log_set:
            model.eval()
            with torch.no_grad():
                ltr, lva = model(Xtr), model(Xva)
                tr_hist[j] = crit(ltr, ytr)
                va_hist[j] = crit(lva, yva)
                if epoch in ckpt_set:
                    thr = (find_mcc_optimal_thresholds(y_tr, torch.sigmoid(ltr).cpu().numpy())
                           if tune_thr else 0.5)
                    ckpt[epoch] = val_mcc_per_label(y_va, torch.sigmoid(lva).cpu().numpy(), thr)
            j += 1
    return tr_hist.cpu().numpy(), va_hist.cpu().numpy(), ckpt


def curve_frame(meta, log_epochs, tr_hist, va_hist):
    n = len(log_epochs)
    return pd.DataFrame({**meta, "epoch": np.repeat(log_epochs, len(LABELS)),
                         "label": np.tile(LABELS, n),
                         "train_loss": tr_hist.reshape(-1), "val_loss": va_hist.reshape(-1)})


def ckpt_frame(meta, ckpt):
    rows = [{**meta, "epoch": ep, "label": lab, "mcc": m}
            for ep, ms in ckpt.items() for lab, m in zip(LABELS, ms)]
    return pd.DataFrame(rows)


# ============================================================================
# Task runner
# ============================================================================

def run(task_id, smoke):
    set_name, embedding, only_loss = TASKS[task_id]
    out_dir = OUT_DIR / "_smoke" if smoke else OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = "__".join(x for x in (set_name, embedding, only_loss) if x)
    f_curves, f_ckpt = out_dir / f"{stem}__curves.csv.gz", out_dir / f"{stem}__ckpt_mcc.csv.gz"
    print(f"task {task_id}: {stem} | device={DEVICE} | smoke={smoke}", flush=True)

    if set_name == "up":
        X, y, splits = up_dataset_and_splits()
        tune_thr, ds_label = False, "UP"
    else:
        ds = set_name.split("_")[1]
        group = GROUP_OF[embedding]
        y, feats, splits = sam_dataset_and_splits(group, ds)
        X = np.hstack([feats[s] for s in SAM_GROUPS[group]["conditions"][embedding]]).astype(np.float32)
        tune_thr, ds_label = True, ds
    print(f"features {X.shape} | label totals {dict(zip(LABELS, y.sum(0).astype(int)))}", flush=True)

    max_epochs, log_epochs, ckpt_epochs, n_iter = MAX_EPOCHS, LOG_EPOCHS, EPOCH_CHECKPOINTS, N_ITERATIONS
    if smoke:
        max_epochs, n_iter = 60, 1
        log_epochs, ckpt_epochs = list(range(1, 61)), [20, 40, 60]

    losses = [only_loss] if only_loss else LOSS_TYPES
    curves = [pd.read_csv(f_curves)] if f_curves.exists() else []
    ckpts = [pd.read_csv(f_ckpt)] if f_ckpt.exists() else []
    done = (set(zip(curves[0]["iteration"], curves[0]["loss_type"])) if curves else set())

    for it in range(n_iter):
        tr, va = splits[it]
        for li, loss_type in enumerate(losses):
            if (it + 1, loss_type) in done:
                continue
            t0 = time.time()
            seed = SEED + 100 * it + LOSS_TYPES.index(loss_type)
            tr_h, va_h, ck = train_with_loss_curves(X[tr], y[tr], X[va], y[va], loss_type, seed,
                                                    tune_thr, max_epochs, log_epochs, ckpt_epochs)
            meta = {"dataset": ds_label, "embedding": embedding, "loss_type": loss_type,
                    "iteration": it + 1}
            curves.append(curve_frame(meta, log_epochs, tr_h, va_h))
            ckpts.append(ckpt_frame(meta, ck))
            pd.concat(curves).to_csv(f_curves, index=False)     # incremental: resumable
            pd.concat(ckpts).to_csv(f_ckpt, index=False)
            tot_va = va_h.sum(1)
            print(f"  iter {it + 1}/{n_iter} {loss_type:18s} {time.time() - t0:6.0f}s | "
                  f"final train {tr_h[-1].sum():.4g} val {tot_va[-1]:.4g} | "
                  f"min val {tot_va.min():.4g} @ epoch {log_epochs[int(tot_va.argmin())]}", flush=True)
    print("DONE", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-id", type=int, default=None,
                    help=f"0..{len(TASKS) - 1} (SLURM_ARRAY_TASK_ID); --list to enumerate")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--smoke", action="store_true",
                    help="60 epochs, 1 iteration, writes to train_val_loss_output/_smoke/")
    a = ap.parse_args()
    if a.list or a.task_id is None:
        for i, t in enumerate(TASKS):
            print(i, *[x for x in t if x])
        return
    torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    run(a.task_id, a.smoke)


if __name__ == "__main__":
    main()
