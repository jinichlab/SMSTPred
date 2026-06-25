#!/usr/bin/env python3
"""
Per-residue Integrated Gradients over ALL no-motif family reps, aggregated to
surface TRENDS in what the classifier keys on when there is no linear motif.

For each no-motif rep we run IG (logit space, mean-MTase baseline) toward its
assigned substrate class, then aggregate across sequences three ways:

  1. AA-identity enrichment  -- of the top-10% IG-positive residues vs background,
     per class. "When there's no motif, does the call lean on Cys / aromatics /
     charged residues?"  (log2 enrichment per amino acid)
  2. Peakedness (Gini)       -- how concentrated each sequence's attribution is.
     Low Gini = genuinely distributed signal (not a motif we merely missed).
  3. High-attribution windows -- the +/-4 window around each sequence's single
     most important residue, pooled per class -> top_windows_{cls}.fasta, ready
     to STREME for short/degenerate sub-motifs STREME's 8-wide search missed.

Inputs:
  no_motif.fasta, no_motif_manifest.csv      (no_motif_prep.py)
  no_motif_prostt5.h5                         (embed_prostt5.py, per-residue ON)
  deployment model + ../corpus.csv (in_domain)
Outputs (./no_motif_attr_output/):
  per_sequence.csv, aa_enrichment.csv, aa_enrichment.png, peakedness.png,
  top_windows_{DNA,Protein,RNA,Other}.fasta

Env: classifier_env (GPU).
"""
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from integrated_gradients import load_ensemble, mean_baseline, integrated_gradients, class_score  # noqa: E402

MODEL_DIR = Path("/home/nroethler/Code/journal_article/repo/2_training/deployment_model")
CORPUS = HERE.parent / "corpus.csv"
RESIDUE_H5 = HERE / "no_motif_prostt5.h5"
MANIFEST = HERE / "no_motif_manifest.csv"
FASTA = HERE / "no_motif.fasta"
OUT = HERE / "no_motif_attr_output"

AA = "ACDEFGHIKLMNPQRSTVWY"
SUBSTRATE_COLORS = {"Other": "#e6550d", "Protein": "#238b45", "RNA": "#cb181d", "DNA": "#6a51a3"}
CLASSES = ["DNA", "Protein", "RNA", "Other"]
TOP_FRAC = 0.10          # top 10% IG-positive residues define "essential"
WIN = 4                  # +/-4 window around the single top residue


def fasta_seqs(path):
    seqs, sid, buf = {}, None, []
    for ln in open(path):
        ln = ln.rstrip()
        if ln.startswith(">"):
            if sid:
                seqs[sid] = "".join(buf)
            sid, buf = ln[1:].split()[0], []
        elif ln:
            buf.append(ln)
    if sid:
        seqs[sid] = "".join(buf)
    return seqs


def gini(x):
    x = np.clip(x, 0, None).astype(np.float64)
    if x.sum() == 0:
        return 0.0
    x = np.sort(x)
    n = len(x)
    return float((2 * np.sum((np.arange(1, n + 1)) * x) / (n * x.sum())) - (n + 1) / n)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    cfg, labels, thresholds, models = load_ensemble(MODEL_DIR, device)
    baseline = mean_baseline(device)
    seqs = fasta_seqs(FASTA)
    man = pd.read_csv(MANIFEST).set_index("Entry")
    in_dom = pd.read_csv(CORPUS, usecols=["Entry", "in_domain"]).drop_duplicates("Entry").set_index("Entry")["in_domain"]
    print(f"device={device} | {len(man)} no-motif reps | labels={labels}")

    # per-class accumulators
    top_cnt = {c: defaultdict(int) for c in CLASSES}     # AA counts in top residues
    bg_cnt = {c: defaultdict(int) for c in CLASSES}      # AA counts overall
    win_out = {c: open(OUT / f"top_windows_{c}.fasta", "w") for c in CLASSES}
    rows = []

    with h5py.File(RESIDUE_H5, "r") as f:
        for n, (entry, r) in enumerate(man.iterrows()):
            cls = r["class"]
            if entry not in f or "residue" not in f[entry]:
                continue
            seq = seqs.get(entry, "")
            E = torch.tensor(np.array(f[entry]["residue"], np.float32), device=device)
            c = labels.index(cls)
            A, _, _, Fz, Fb, comp = integrated_gradients(models, E, c, baseline, 64, device, "logit")
            L = len(A)

            # 1) AA enrichment: top-10% IG-positive residues
            k = max(1, int(TOP_FRAC * L))
            top_pos = np.argsort(A)[-k:]
            for i in range(L):
                aa = seq[i] if i < len(seq) else "X"
                if aa in AA:
                    bg_cnt[cls][aa] += 1
            for i in top_pos:
                aa = seq[i] if i < len(seq) else "X"
                if aa in AA:
                    top_cnt[cls][aa] += 1

            # 3) window around the single most-essential residue
            p = int(np.argmax(A))
            w = seq[max(0, p - WIN): p + WIN + 1]
            if len(w) >= 5:
                win_out[cls].write(f">{entry}_p{p+1}\n{w}\n")

            g = gini(A)
            top5 = np.argsort(A)[-5:][::-1]
            rows.append({"Entry": entry, "class": cls, "length": L,
                         "in_domain": int(in_dom.get(entry, -1)),
                         "gini": round(g, 4), "completeness": round(comp[0] / comp[1], 3) if abs(comp[1]) > 1e-9 else np.nan,
                         "top_residues": ";".join(f"{(seq[i] if i < len(seq) else '?')}{i+1}" for i in top5)})
            if (n + 1) % 500 == 0:
                print(f"  {n+1}/{len(man)} ...", flush=True)

    for fh in win_out.values():
        fh.close()
    per = pd.DataFrame(rows)
    per.to_csv(OUT / "per_sequence.csv", index=False)
    print(f"completeness (should be ~1.0): {per.completeness.median():.3f} median")

    # ---- AA enrichment table + figure ----
    enr_rows = []
    for c in CLASSES:
        tt = sum(top_cnt[c].values()); bt = sum(bg_cnt[c].values())
        for aa in AA:
            tf = top_cnt[c][aa] / tt if tt else 0
            bf = bg_cnt[c][aa] / bt if bt else 0
            enr_rows.append({"class": c, "aa": aa, "top_freq": tf, "bg_freq": bf,
                             "log2_enrich": np.log2((tf + 1e-6) / (bf + 1e-6))})
    enr = pd.DataFrame(enr_rows)
    enr.to_csv(OUT / "aa_enrichment.csv", index=False)

    fig, axes = plt.subplots(1, 4, figsize=(20, 4.4), sharey=True)
    for ax, c in zip(axes, CLASSES):
        sub = enr[enr["class"] == c].set_index("aa").loc[list(AA)]
        cols = [SUBSTRATE_COLORS[c] if v > 0 else "#9e9e9e" for v in sub.log2_enrich]
        ax.bar(range(len(AA)), sub.log2_enrich, color=cols, edgecolor="black", linewidth=0.4)
        ax.set_xticks(range(len(AA))); ax.set_xticklabels(list(AA), fontsize=8)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(f"{c}  (n={ (per['class']==c).sum() })", fontweight="bold")
        ax.grid(True, axis="y", alpha=0.3)
    axes[0].set_ylabel("log2 enrichment in top-10% IG residues")
    fig.suptitle("No-motif reps — amino-acid enrichment among the residues the classifier "
                 "relies on most (IG, logit space)", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(OUT / "aa_enrichment.png", dpi=150)
    plt.close(fig)

    # ---- peakedness (Gini) per class ----
    fig, ax = plt.subplots(figsize=(7, 4.4))
    data = [per[per["class"] == c]["gini"].values for c in CLASSES]
    parts = ax.violinplot(data, showmedians=True)
    ax.set_xticks(range(1, len(CLASSES) + 1)); ax.set_xticklabels(CLASSES)
    ax.set_ylabel("Gini of |IG| (0 = distributed, 1 = focal)")
    ax.set_title("How concentrated is the no-motif attribution signal?", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(OUT / "peakedness.png", dpi=150); plt.close(fig)

    print("\n=== AA enrichment: top-3 enriched residues per class (top-10% IG) ===")
    for c in CLASSES:
        s = enr[enr["class"] == c].sort_values("log2_enrich", ascending=False).head(3)
        print(f"  {c:8} " + ", ".join(f"{r.aa}(+{r.log2_enrich:.2f})" for _, r in s.iterrows())
              + f"   | median Gini={per[per['class']==c]['gini'].median():.3f}")
    print(f"\n[done] -> {OUT}")


if __name__ == "__main__":
    main()
