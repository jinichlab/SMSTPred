#!/usr/bin/env python3
"""
SAM-trained vs UP-trained ProstT5 MLP, evaluated on the SAME held-out SAM test set.

Question: for SAM-MTase substrate prediction, does training on the broad 201k UP
enzyme set match focused SAM training? Same ProstT5 MLP architecture; only the
training data differs.

Design (MVP -- exact-entry holdout; homolog exclusion deferred)
--------------------------------------------------------------
SAM is clustered once (greedy cosine on ProtT5 @ 0.6, in-script like the deployment).
Over N=16 iterations, whole clusters are shuffled into a fresh 90/10 train/test split
(cluster-disjoint within SAM; the split is retried until every label, including the
108-protein DNA class, has a positive in the test set). Per iteration:
  TEST     = the held-out 10% of SAM clusters (cluster-disjoint within SAM).
  SAM model trains on the 90% train pool (cluster-disjoint).
  UP  model trains on ALL UP proteins MINUS the exact TEST entry IDs.
  Both predict the SAME held-out SAM test set.
Metrics are mean +/- SD over the 16 repeated splits. This mirrors Panel A's
16-iteration repeated-split design, so the two Figure-2 panels are consistent.

  *** CAVEAT (known, deferred): the UP training set still contains the non-test SAM
  proteins AND any UP homologs of the test proteins -- only the exact test entries are
  removed. So the UP model has a mild homology advantage over the cluster-disjoint SAM
  model. A later version will also drop UP proteins within 0.6 cosine of the test set. ***

Both models use their own deployment config:
  SAM:  300 epochs, loss {Other: sqrt_weighted_bce, Protein/RNA/DNA: bce}
  UP : 3500 epochs, loss {Other/Protein: weighted_bce, RNA/DNA: sqrt_weighted_bce}
Per-label decision thresholds are MCC-optimal, tuned on each model's own TRAIN
predictions (matching the deployment), so MCC/F1 are at the deployment operating point.

Outputs (./sam_vs_up_prostt5_comparison_output/)
  sam_vs_up_prostt5_comparison_summary.csv    model x label x metric, mean/std over iters
  sam_vs_up_prostt5_comparison_detailed.csv    per-iteration rows
  sam_vs_up_prostt5_comparison_mcc.png         per-label MCC + AUPRC, SAM vs UP

Env: classifier_env (GPU).
"""
import os
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import (matthews_corrcoef, f1_score, accuracy_score,
                             precision_recall_curve, auc, roc_auc_score)
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

random.seed(42); np.random.seed(42); torch.manual_seed(42)

FILES = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[1] / "data"))
SAM_CSV = FILES / "labels_SAM_full.csv"
UP_CSV = FILES / "labels_UP.csv"
PROSTT5_SAM = FILES / "embeddings_SAM_ProstT5.h5"     # SAM features (mean group)
PROTT5_SAM = FILES / "embeddings_SAM_ProtT5.h5"              # SAM clustering space (flat)
PROSTT5_UP = FILES / "embeddings_UP_ProstT5.h5"       # UP features (mean group)

OUT = Path(__file__).resolve().parent / "sam_vs_up_prostt5_comparison_output"
LABELS = ["Other", "Protein", "RNA", "DNA"]
CLUSTER_THRESHOLD = 0.6
N_ITER = 16              # repeated cluster-disjoint 90/10 splits (matches Panel A)
TRAIN_RATIO = 0.9
LR = 0.001
HIDDEN = (512, 256)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

SAM_CFG = {"epochs": 300,
           "loss": {0: "sqrt_weighted_bce", 1: "bce", 2: "bce", 3: "bce"}}
UP_CFG = {"epochs": 3500,
          "loss": {0: "weighted_bce", 1: "weighted_bce", 2: "sqrt_weighted_bce", 3: "sqrt_weighted_bce"}}


# ---------------------------------------------------------------- data loading
def load_mean_group(path):
    with h5py.File(path, "r") as f:
        return {k: np.asarray(f[k]["mean"], np.float32) for k in f.keys()}


def load_flat(path):
    with h5py.File(path, "r") as f:
        return {k: np.asarray(f[k], np.float32) for k in f.keys()}


def load_sam():
    prostt5 = load_mean_group(PROSTT5_SAM)
    prott5 = load_flat(PROTT5_SAM)
    df = pd.read_csv(SAM_CSV, low_memory=False).drop_duplicates("Entry")
    df = df[df["Entry"].apply(lambda e: e in prostt5 and e in prott5)].reset_index(drop=True)
    entries = df["Entry"].values
    Xp = np.vstack([prostt5[e] for e in entries]).astype(np.float32)
    Xc = np.vstack([prott5[e] for e in entries]).astype(np.float32)
    y = df[LABELS].values.astype(np.float32)
    return entries, Xp, Xc, y


def load_up():
    prostt5 = load_mean_group(PROSTT5_UP)
    df = pd.read_csv(UP_CSV, low_memory=False).drop_duplicates("Entry")
    df = df[df["Entry"].apply(lambda e: e in prostt5)].reset_index(drop=True)
    entries = df["Entry"].values
    Xp = np.vstack([prostt5[e] for e in entries]).astype(np.float32)
    y = df[LABELS].values.astype(np.float32)
    return entries, Xp, y


# ------------------------------------------------------- clustering + folds
def greedy_clustering(sim, threshold):
    n = len(sim)
    order = np.argsort(-np.linalg.norm(sim, axis=1))
    clusters, labels, cid = {}, -np.ones(n, dtype=int), 0
    for i in tqdm(order, desc="cluster SAM@0.6", leave=False):
        if labels[i] != -1:
            continue
        for rep in clusters:
            if sim[i, clusters[rep][0]] >= threshold:
                labels[i] = rep; clusters[rep].append(i); break
        else:
            clusters[cid] = [i]; labels[i] = cid; cid += 1
    return clusters


def split_clusters_train_test(clusters, y, train_ratio=TRAIN_RATIO, max_attempts=1000):
    """Shuffle whole clusters into a train/test split; retry until every label has a
    positive in test (guards the rare DNA class -> MCC always defined). Mirrors
    Panel A's split_clusters_train_test. Relies on the module RNG state, so calling
    it repeatedly in a loop yields distinct-but-reproducible splits."""
    cluster_ids = list(clusters.keys())
    n_labels = y.shape[1]
    for _ in range(max_attempts):
        random.shuffle(cluster_ids)
        k = int(train_ratio * len(cluster_ids))
        train_idx = [i for c in cluster_ids[:k] for i in clusters[c]]
        test_idx = [i for c in cluster_ids[k:] for i in clusters[c]]
        if all(y[test_idx][:, j].sum() > 0 for j in range(n_labels)):
            return np.array(train_idx), np.array(test_idx)
    raise ValueError("Could not find a valid cluster split with all labels in test.")


# ----------------------------------------------------------- model + loss
class MultilabelNN(nn.Module):
    def __init__(self, in_dim, out_dim, hidden=HIDDEN):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hidden[0]); self.fc2 = nn.Linear(hidden[0], hidden[1])
        self.fc3 = nn.Linear(hidden[1], out_dim); self.relu = nn.ReLU()

    def forward(self, x):
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))


def pos_weight(y):
    c = y.sum(0).clip(min=1)
    return torch.tensor((len(y) - c) / c, dtype=torch.float32)


class PerLabelMixedLoss(nn.Module):
    def __init__(self, pw, loss_map, gamma=2.0):
        super().__init__(); self.pw = pw; self.lm = loss_map; self.gamma = gamma

    def forward(self, logits, targets):
        total = 0.0
        for i in range(logits.shape[1]):
            lt = self.lm.get(i, "bce"); li, ti = logits[:, i], targets[:, i]
            if lt == "focal":
                bce = nn.functional.binary_cross_entropy_with_logits(li, ti, reduction="none")
                p = torch.sigmoid(li); p_t = ti * p + (1 - ti) * (1 - p)
                w = ti * self.pw[i] + (1 - ti)
                total += (((1 - p_t) ** self.gamma) * bce * w).mean()
            elif lt == "sqrt_weighted_bce":
                total += nn.functional.binary_cross_entropy_with_logits(li, ti, pos_weight=torch.sqrt(self.pw[i:i + 1]))
            elif lt == "weighted_bce":
                total += nn.functional.binary_cross_entropy_with_logits(li, ti, pos_weight=self.pw[i:i + 1])
            else:
                total += nn.functional.binary_cross_entropy_with_logits(li, ti)
        return total


def train_mlp(X, y, cfg, seed):
    torch.manual_seed(seed)
    model = MultilabelNN(X.shape[1], y.shape[1]).to(DEVICE)
    crit = PerLabelMixedLoss(pos_weight(y).to(DEVICE), cfg["loss"])
    opt = optim.Adam(model.parameters(), lr=LR)
    Xt = torch.tensor(X, dtype=torch.float32, device=DEVICE)
    yt = torch.tensor(y, dtype=torch.float32, device=DEVICE)
    n = X.shape[0]; k = max(1, int(0.9 * n))
    for _ in range(cfg["epochs"]):
        idx = torch.randperm(n, device=DEVICE)[:k]
        model.train(); opt.zero_grad()
        crit(model(Xt[idx]), yt[idx]).backward(); opt.step()
    return model


@torch.no_grad()
def predict(model, X):
    model.eval()
    out = []
    for i in range(0, len(X), 8192):
        xb = torch.tensor(X[i:i + 8192], dtype=torch.float32, device=DEVICE)
        out.append(torch.sigmoid(model(xb)).cpu().numpy())
    return np.vstack(out)


def optimal_thresholds(y, probs):
    thr = np.full(len(LABELS), 0.5)
    for i in range(len(LABELS)):
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


def metrics(y, probs, thr):
    out = {}
    for i, lab in enumerate(LABELS):
        pred = (probs[:, i] >= thr[i]).astype(int)
        rec = {"mcc": matthews_corrcoef(y[:, i], pred) if len(np.unique(y[:, i])) > 1 else np.nan,
               "f1": f1_score(y[:, i], pred, zero_division=0),
               "accuracy": accuracy_score(y[:, i], pred)}
        if len(np.unique(y[:, i])) > 1:
            pr, rc, _ = precision_recall_curve(y[:, i], probs[:, i])
            rec["auprc"] = auc(rc, pr); rec["auroc"] = roc_auc_score(y[:, i], probs[:, i])
        else:
            rec["auprc"] = rec["auroc"] = np.nan
        out[lab] = rec
    return out


# ---------------------------------------------------------------------- main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("=" * 70 + f"\nSAM vs UP ProstT5 MLP  (device={DEVICE})\n" + "=" * 70)

    sam_e, Xsam_p, Xsam_c, ysam = load_sam()
    up_e, Xup_p, yup = load_up()
    print(f"SAM {len(sam_e)} | UP {len(up_e)} | "
          f"SAM labels {dict(zip(LABELS, ysam.sum(0).astype(int)))}")

    clusters = greedy_clustering(cosine_similarity(Xsam_c), CLUSTER_THRESHOLD)
    print(f"{len(clusters)} SAM clusters -> {N_ITER} repeated {int(TRAIN_RATIO*100)}/"
          f"{int((1-TRAIN_RATIO)*100)} cluster-disjoint splits")
    up_idx = {e: i for i, e in enumerate(up_e)}

    rows = []
    for it in range(N_ITER):
        sam_tr, te = split_clusters_train_test(clusters, ysam)
        test_set = set(sam_e[te])
        up_tr = np.array([i for e, i in up_idx.items() if e not in test_set])
        print(f"\niter {it}: test {len(te)} | SAM train {len(sam_tr)} | "
              f"UP train {len(up_tr)} (removed {len(up_e) - len(up_tr)} test entries) | "
              f"test labels {dict(zip(LABELS, ysam[te].sum(0).astype(int)))}")

        for name, (Xtr, ytr, cfg, seed) in {
            "SAM": (Xsam_p[sam_tr], ysam[sam_tr], SAM_CFG, it),
            "UP":  (Xup_p[up_tr],   yup[up_tr],   UP_CFG, 100 + it),
        }.items():
            model = train_mlp(Xtr, ytr, cfg, seed)
            thr = optimal_thresholds(ytr, predict(model, Xtr))   # tuned on TRAIN
            mt = metrics(ysam[te], predict(model, Xsam_p[te]), thr)
            for lab in LABELS:
                rows.append({"model": name, "iteration": it, "label": lab, **mt[lab]})
            mcc_mean = np.nanmean([mt[l]["mcc"] for l in LABELS])
            print(f"  {name:3} thr={np.round(thr,3)}  meanMCC={mcc_mean:.4f}  "
                  + " ".join(f"{l}:{mt[l]['mcc']:.3f}" for l in LABELS))

    det = pd.DataFrame(rows)
    det.to_csv(OUT / "sam_vs_up_prostt5_comparison_detailed.csv", index=False)

    METR = ["mcc", "auprc", "auroc", "f1", "accuracy"]
    summ = (det.groupby(["model", "label"])[METR].agg(["mean", "std"]).reset_index())
    summ.columns = ["model", "label"] + [f"{m}_{s}" for m in METR for s in ("mean", "std")]
    # overall = mean across labels, per model/iteration, then mean/std over iterations
    ov = det.groupby(["model", "iteration"])[METR].mean().reset_index()
    ovs = ov.groupby("model")[METR].agg(["mean", "std"]).reset_index()
    ovs.columns = ["model"] + [f"{m}_{s}" for m in METR for s in ("mean", "std")]
    ovs.insert(1, "label", "Overall")
    summ = pd.concat([summ, ovs], ignore_index=True)
    summ.to_csv(OUT / "sam_vs_up_prostt5_comparison_summary.csv", index=False)

    print("\n=== SAM vs UP  (mean +/- SD over iterations, tuned-threshold) ===")
    for lab in LABELS + ["Overall"]:
        r = {m: summ[(summ.label == lab) & (summ.model == m)] for m in ["SAM", "UP"]}
        s_mcc = r["SAM"].mcc_mean.values[0]; u_mcc = r["UP"].mcc_mean.values[0]
        s_ap = r["SAM"].auprc_mean.values[0]; u_ap = r["UP"].auprc_mean.values[0]
        print(f"  {lab:8} MCC  SAM={s_mcc:.4f}  UP={u_mcc:.4f}  (Δ{u_mcc - s_mcc:+.4f}) | "
              f"AUPRC SAM={s_ap:.4f} UP={u_ap:.4f}")

    # ---- figure: per-label MCC and AUPRC, SAM vs UP ----
    order = LABELS + ["Overall"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, metric, lbl in zip(axes, ["mcc", "auprc"], ["MCC (tuned threshold)", "AUPRC"]):
        x = np.arange(len(order)); w = 0.38
        for k, (mdl, col) in enumerate([("SAM", "#2c7fb8"), ("UP", "#e6550d")]):
            means = [summ[(summ.label == l) & (summ.model == mdl)][f"{metric}_mean"].values[0] for l in order]
            stds = [summ[(summ.label == l) & (summ.model == mdl)][f"{metric}_std"].values[0] for l in order]
            ax.bar(x + (k - 0.5) * w, means, w, yerr=stds, capsize=3, label=f"{mdl}-trained",
                   color=col, edgecolor="black", linewidth=0.5, error_kw=dict(lw=1, ecolor="#444"))
        ax.set_xticks(x); ax.set_xticklabels(order); ax.set_ylabel(lbl)
        ax.set_ylim(0, 1.0); ax.grid(True, axis="y", alpha=0.3); ax.legend()
        ax.set_title(lbl, fontweight="bold")
    fig.suptitle("SAM-trained vs UP-trained ProstT5 MLP — evaluated on held-out SAM "
                 f"({N_ITER} cluster-disjoint splits, mean±SD)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "sam_vs_up_prostt5_comparison_mcc.png", dpi=150)
    print(f"\nwrote summary/detailed CSVs + figure -> {OUT}")


if __name__ == "__main__":
    main()
