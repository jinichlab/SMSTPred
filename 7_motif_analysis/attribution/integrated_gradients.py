#!/usr/bin/env python3
"""
Per-residue Integrated Gradients for the ProstT5 -> MLP substrate classifier.

WHY THIS, NOT PLAIN SALIENCY
----------------------------
The classifier eats a MEAN-POOLED ProstT5 embedding:  z = (1/L) * sum_i e_i ,
logit = f(z). Because pooling is linear, d(logit)/d(e_i) is IDENTICAL for every
residue -> naive gradient saliency is uniform and useless. We instead attribute
through each residue's OWN embedding e_i, which differ.

METHOD
------
Target function F_c(z) = ensemble-mean probability of class c
                       = mean_m sigmoid( member_m(z)[c] )      (the deployed score).

Integrated Gradients along the pooled-embedding path from a baseline b to z:
    G = (1/S) * sum_{k} grad_z F_c( b + a_k (z - b) ) ,  a_k = (k-0.5)/S  (midpoint)
Per-residue attribution (zero or constant baseline b_res, pooling to b):
    A_i = (1/L) * (e_i - b_res) . G
Completeness:  sum_i A_i = (z - b).G = F_c(z) - F_c(b)   (checked & reported).

A_i > 0  => residue i pushes the class-c probability UP relative to baseline.

BASELINE
--------
  --baseline mean (default): b = mean of all training "mean" embeddings, i.e. a
     "typical MTase". Reads as "what makes THIS enzyme a <class>-MTase rather than
     a generic one" -- the substrate-discriminative signal we care about.
  --baseline zero: b = 0 (standard IG); reads as absolute contribution. Off-manifold.

Also emits raw Gradient x Input (single gradient at z) as a cheap companion column.

INPUT  : an H5 with per-residue ProstT5 embeddings ( /<Entry>/residue (L,1024),
         /<Entry>/mean (1024,) ), as produced by embed_prostt5.py (per-residue ON).
OUTPUT : per-entry residue-attribution CSV + PNG track, and a run summary CSV.

Env: classifier_env.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path("/home/nroethler/Code/journal_article/repo")
DEFAULT_MODEL = REPO / "2_training" / "deployment_model"
AD_REF_H5 = Path("/home/nroethler/Code/files/prostt5_SAM_swiss.h5")  # for "mean" baseline
SUBSTRATE_COLORS = {"Other": "#e6550d", "Protein": "#238b45",
                    "RNA": "#cb181d", "DNA": "#6a51a3"}


class MultilabelNN(nn.Module):
    """Must match train_deployment_ensemble.py exactly (1024->512->256->4, logits)."""
    def __init__(self, in_dim, out_dim, hidden):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hidden[0])
        self.fc2 = nn.Linear(hidden[0], hidden[1])
        self.fc3 = nn.Linear(hidden[1], out_dim)
        self.relu = nn.ReLU()

    def forward(self, x):
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))


def load_ensemble(model_dir, device):
    cfg = json.load(open(model_dir / "model_config.json"))
    thr = json.load(open(model_dir / "thresholds.json"))["thresholds"]
    labels = cfg["labels"]
    models = []
    for s in range(cfg["n_ensemble"]):
        m = MultilabelNN(cfg["input_dim"], len(labels), cfg["hidden"]).to(device)
        m.load_state_dict(torch.load(model_dir / f"ensemble_member_{s}.pt", map_location=device))
        m.eval()
        for p in m.parameters():
            p.requires_grad_(False)
        models.append(m)
    return cfg, labels, np.array([thr[l] for l in labels], np.float32), models


def class_score(models, z, c, space="prob"):
    """Ensemble target for class index c. z: (..., 1024) -> (...,).

    space="prob"  : mean_m sigmoid(logit) -- the deployed score, but saturates at
                    high confidence (gradient -> 0), and the mean baseline is
                    degenerate for the majority class (RNA).
    space="logit" : mean_m logit -- unsaturated, informative gradients for every
                    class incl. the majority; monotone proxy for the deployed score.
    """
    per = torch.stack([m(z)[..., c] for m in models])  # (M, ...) raw logits
    if space == "prob":
        per = torch.sigmoid(per)
    return per.mean(0)


def mean_baseline(device):
    with h5py.File(AD_REF_H5, "r") as f:
        vals = [np.array(f[k]["mean"], np.float32) for k in f.keys()]
    return torch.tensor(np.mean(vals, axis=0), dtype=torch.float32, device=device)


def integrated_gradients(models, E, c, baseline, steps, device, space="prob"):
    """
    E: (L, 1024) residue embeddings. Returns (A_ig, A_gxi, z, completeness).
    A_*: (L,) per-residue attribution to class index c.
    """
    L = E.shape[0]
    z = E.mean(0)                                   # (1024,)
    b = baseline                                    # (1024,)
    # accumulated IG gradient along the pooled path (midpoint Riemann)
    alphas = (torch.arange(steps, device=device).float() + 0.5) / steps
    path = b[None] + alphas[:, None] * (z - b)[None]    # (S, 1024)
    path.requires_grad_(True)
    F = class_score(models, path, c, space).sum()
    g_path, = torch.autograd.grad(F, path)              # (S, 1024)
    G = g_path.mean(0)                                  # (1024,)  averaged gradient
    A_ig = ((E - b) @ G) / L                            # (L,)

    # raw gradient x input (single gradient at z)
    z1 = z.clone().requires_grad_(True)
    f1 = class_score(models, z1[None], c, space).sum()
    g1, = torch.autograd.grad(f1, z1)
    A_gxi = ((E - b) @ g1) / L

    with torch.no_grad():
        Fz = float(class_score(models, z[None], c, space))
        Fb = float(class_score(models, b[None], c, space))
    completeness = (float(A_ig.sum()), Fz - Fb)         # should match
    return (A_ig.detach().cpu().numpy(), A_gxi.detach().cpu().numpy(),
            z.detach().cpu().numpy(), Fz, Fb, completeness)


def plot_track(entry, seq, A, label, color, score, out_png, topk=8):
    L = len(A)
    x = np.arange(1, L + 1)
    fig, ax = plt.subplots(figsize=(max(8, L / 25), 3.2))
    pos = A.copy(); pos[pos < 0] = 0
    neg = A.copy(); neg[neg > 0] = 0
    ax.bar(x, pos, width=1.0, color=color, linewidth=0)
    ax.bar(x, neg, width=1.0, color="#9e9e9e", linewidth=0)
    ax.axhline(0, color="k", lw=0.6)
    # annotate top-|A| positions with residue letter
    top = np.argsort(-np.abs(A))[:topk]
    for i in sorted(top):
        aa = seq[i] if seq and i < len(seq) else ""
        ax.annotate(f"{aa}{i+1}", (i + 1, A[i]),
                    textcoords="offset points", xytext=(0, 3 if A[i] >= 0 else -9),
                    ha="center", fontsize=7,
                    color=color if A[i] >= 0 else "#555555")
    ax.set_xlim(0.5, L + 0.5)
    ax.set_xlabel("residue position")
    ax.set_ylabel(f"IG attribution -> {label}")
    ax.set_title(f"{entry}   {label} score={score:.3f}   "
                 f"(+ = pushes toward {label})", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def fasta_seqs(path):
    seqs, sid, buf = {}, None, []
    for ln in open(path):
        ln = ln.rstrip()
        if not ln:
            continue
        if ln.startswith(">"):
            if sid:
                seqs[sid] = "".join(buf)
            sid, buf = ln[1:].split()[0], []
        else:
            buf.append(ln)
    if sid:
        seqs[sid] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--residue-h5", required=True, type=Path,
                    help="H5 with /<Entry>/residue + /mean (embed_prostt5.py, per-residue ON)")
    ap.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL)
    ap.add_argument("--fasta", type=Path, default=None,
                    help="optional FASTA to annotate residue letters on plots")
    ap.add_argument("--entries", default=None, help="comma-separated subset (default: all in H5)")
    ap.add_argument("--label", default="auto",
                    help="auto (top called/scored), all, or one of Other/Protein/RNA/DNA")
    ap.add_argument("--baseline", choices=["mean", "zero"], default="mean")
    ap.add_argument("--space", choices=["logit", "prob"], default="logit",
                    help="attribute logit (unsaturated, works for majority class) "
                         "or prob (deployed score, saturates at high confidence)")
    ap.add_argument("--steps", type=int, default=64)
    ap.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent / "ig_output")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    cfg, labels, thresholds, models = load_ensemble(args.model_dir, device)
    print(f"device={device} | labels={labels} | baseline={args.baseline} | steps={args.steps}")

    baseline = (mean_baseline(device) if args.baseline == "mean"
                else torch.zeros(cfg["input_dim"], device=device))
    seqs = fasta_seqs(args.fasta) if args.fasta else {}

    with h5py.File(args.residue_h5, "r") as f:
        all_entries = list(f.keys())
        want = ([e.strip() for e in args.entries.split(",")] if args.entries else all_entries)
        summary = []
        for e in want:
            if e not in f:
                print(f"  [warn] {e} not in H5, skipped"); continue
            if "residue" not in f[e]:
                print(f"  [warn] {e} has no /residue (re-embed with per-residue ON), skipped"); continue
            E = torch.tensor(np.array(f[e]["residue"], np.float32), device=device)
            z = E.mean(0)
            with torch.no_grad():
                scores = {lab: float(class_score(models, z[None], i, "prob"))
                          for i, lab in enumerate(labels)}
            called = [lab for i, lab in enumerate(labels) if scores[lab] >= thresholds[i]]

            if args.label == "all":
                targets = labels
            elif args.label == "auto":
                targets = [max(called, key=lambda l: scores[l])] if called \
                    else [max(labels, key=lambda l: scores[l])]
            else:
                targets = [args.label]

            seq = seqs.get(e, "")
            for lab in targets:
                c = labels.index(lab)
                A_ig, A_gxi, zc, Fz, Fb, comp = integrated_gradients(
                    models, E, c, baseline, args.steps, device, args.space)
                tag = f"{e}__{lab}"
                pd.DataFrame({
                    "position": np.arange(1, len(A_ig) + 1),
                    "residue": list(seq) if len(seq) == len(A_ig) else [""] * len(A_ig),
                    "ig": A_ig, "grad_x_input": A_gxi,
                    "ig_norm": A_ig / (np.abs(A_ig).max() + 1e-12),
                }).to_csv(args.out_dir / f"{tag}_residue_attr.csv", index=False)
                if seq or True:
                    plot_track(e, seq, A_ig, lab, SUBSTRATE_COLORS.get(lab, "#333333"),
                               scores[lab], args.out_dir / f"{tag}_track.png")
                top = np.argsort(-A_ig)[:6]
                top_str = ";".join(f"{(seq[i] if i < len(seq) else '?')}{i+1}:{A_ig[i]:+.3f}"
                                   for i in top)
                ratio = comp[0] / comp[1] if abs(comp[1]) > 1e-9 else float("nan")
                summary.append({"Entry": e, "label": lab, "score": round(scores[lab], 4),
                                "called": ",".join(called) or "(none)",
                                "F_baseline": round(Fb, 4),
                                "sum_ig": round(comp[0], 4), "Fz_minus_Fb": round(comp[1], 4),
                                "completeness_ratio": round(ratio, 3),
                                "top_positive_residues": top_str})
                print(f"  {tag:24} score={scores[lab]:.3f}  "
                      f"completeness={ratio:.3f}  top: {top_str}")

    pd.DataFrame(summary).to_csv(args.out_dir / "ig_summary.csv", index=False)
    print(f"[done] {len(summary)} attributions -> {args.out_dir}")


if __name__ == "__main__":
    main()
