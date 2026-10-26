#!/usr/bin/env python3
"""
Statistical comparison of classifier MCC across cluster-disjoint CV iterations.

Design: complete paired matrix (N iterations x k classifiers), one MCC per
(iteration, classifier) on the SAME held-out fold -> Demsar (2006) procedure:

  1. Friedman omnibus test (nonparametric repeated-measures on within-fold ranks)
  2. Post-hoc Nemenyi (all-pairs) + critical-difference (CD) diagram
  3. Post-hoc Bonferroni-Dunn vs a control method (MLP_ProstT5)
  4. Headline pairwise Wilcoxon signed-rank (+ Holm correction)

Run per (split_mode, dataset). Cluster-disjoint = the honest evaluation.
"""
import os, itertools
import numpy as np
import pandas as pd
from scipy import stats
import scikit_posthocs as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
# BAKEOFF_OUT (optional) points at another bake-off output folder, e.g. a re-run, without touching the committed one
BASE = os.path.join(HERE, os.environ.get("BAKEOFF_OUT", "sam_classifier_comparison_output"))
SRC = os.path.join(BASE, "sam_classifier_comparison_detailed.csv")
OUT = os.path.join(BASE, "stats")
os.makedirs(OUT, exist_ok=True)
CONTROL = "MLP_ProstT5"
ALPHA = 0.05


def build_matrix(df, split_mode, dataset):
    sub = df[(df.split_mode == split_mode) & (df.dataset == dataset)]
    piv = sub.pivot_table(index="iteration", columns="classifier", values="mcc")
    # order columns by mean MCC (best first)
    piv = piv[piv.mean().sort_values(ascending=False).index]
    return piv  # rows = iterations (blocks), cols = methods


def avg_ranks(piv):
    # rank within each fold: best MCC -> rank 1
    ranks = piv.rank(axis=1, ascending=False)
    return ranks.mean(axis=0).sort_values()  # low rank = better


def nemenyi_cd(k, N, alpha=0.05):
    q = stats.studentized_range.ppf(1 - alpha, k, np.inf) / np.sqrt(2)
    return q * np.sqrt(k * (k + 1) / (6.0 * N))


def bonferroni_dunn_cd(k, N, alpha=0.05):
    # control vs (k-1) others; two-tailed z with Bonferroni over k-1 comparisons
    z = stats.norm.ppf(1 - alpha / (2 * (k - 1)))
    return z * np.sqrt(k * (k + 1) / (6.0 * N))


def cd_diagram(ranks, cd, title, path):
    methods = ranks.index.tolist()
    vals = ranks.values
    k = len(methods)
    lo, hi = int(np.floor(vals.min())), int(np.ceil(vals.max()))
    fig, ax = plt.subplots(figsize=(10, 0.5 * k + 2))
    ax.set_xlim(lo - 0.3, hi + 0.3); ax.set_ylim(0, k + 2)
    ax.invert_xaxis()  # rank 1 (best) on the right-> put best on left instead
    # axis line
    yaxis = k + 1
    ax.hlines(yaxis, lo, hi, color="k")
    for x in range(lo, hi + 1):
        ax.vlines(x, yaxis, yaxis + 0.15, color="k")
        ax.text(x, yaxis + 0.25, str(x), ha="center", va="bottom", fontsize=9)
    # method labels stacked down the sides
    for i, (m, v) in enumerate(zip(methods, vals)):
        y = yaxis - 1 - i
        ax.plot([v, v], [yaxis, y], color="#444", lw=1)
        ax.plot([v, lo if v < (lo + hi) / 2 else hi], [y, y], color="#444", lw=1)
        side = lo if v < (lo + hi) / 2 else hi
        ha = "right" if side == lo else "left"
        ax.text(side + (-0.05 if ha == "right" else 0.05), y,
                f"{m}  ({v:.2f})", ha=ha, va="center", fontsize=9)
    # CD bar
    ax.hlines(yaxis + 1.1, hi, hi - cd, color="crimson", lw=3)
    ax.text(hi - cd / 2, yaxis + 1.25, f"CD = {cd:.2f}", ha="center",
            color="crimson", fontsize=10)
    # cliques: connect methods that are NOT significantly different (within CD)
    y0 = 0.6
    used = np.zeros(k)
    groups = []
    for i in range(k):
        j = i
        while j + 1 < k and (vals[j + 1] - vals[i]) <= cd:
            j += 1
        if j > i:
            groups.append((i, j))
    # merge nested
    filtered = [g for idx, g in enumerate(groups)
                if not any(g2[0] <= g[0] and g[1] <= g2[1] and g2 != g
                           for g2 in groups)]
    for gi, (a, b) in enumerate(filtered):
        yy = y0 + 0.22 * gi
        ax.hlines(yy, vals[a] - 0.03, vals[b] + 0.03, color="steelblue", lw=4)
    ax.set_title(title, fontsize=12)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def analyze(df, split_mode, dataset, log):
    piv = build_matrix(df, split_mode, dataset)
    N, k = piv.shape
    tag = f"{split_mode}/{dataset}"
    log.append(f"\n{'='*70}\n{tag}   (N={N} iterations, k={k} methods)\n{'='*70}")

    # 1. Friedman
    chi2, p = stats.friedmanchisquare(*[piv[c].values for c in piv.columns])
    log.append(f"Friedman omnibus:  chi2={chi2:.3f}, df={k-1}, p={p:.3e}"
               f"   -> {'reject H0 (methods differ)' if p < ALPHA else 'no evidence of difference'}")

    ranks = avg_ranks(piv)
    log.append("\nAverage ranks (1 = best):")
    for m, r in ranks.items():
        log.append(f"   {r:5.2f}  {m:22s}  MCC={piv[m].mean():.3f}+-{piv[m].std():.3f}")

    cd_n = nemenyi_cd(k, N, ALPHA)
    cd_bd = bonferroni_dunn_cd(k, N, ALPHA)
    log.append(f"\nNemenyi critical difference (all-pairs, a=0.05): CD = {cd_n:.3f} rank units")
    log.append(f"Bonferroni-Dunn CD vs control ({CONTROL}):          CD = {cd_bd:.3f} rank units")

    # 2. Nemenyi all-pairs p-value matrix
    nem = sp.posthoc_nemenyi_friedman(piv.values)
    nem.index = piv.columns; nem.columns = piv.columns
    nem.to_csv(os.path.join(OUT, f"nemenyi_{split_mode}_{dataset}.csv"))
    sig = [(a, b, nem.loc[a, b]) for a, b in itertools.combinations(piv.columns, 2)
           if nem.loc[a, b] < ALPHA]
    log.append(f"\nNemenyi significant pairs (p<{ALPHA}): {len(sig)} of {k*(k-1)//2}")
    for a, b, pv in sorted(sig, key=lambda x: x[2]):
        log.append(f"   p={pv:.3f}  {a}  vs  {b}")

    # 3. Bonferroni-Dunn vs control (rank-difference test)
    log.append(f"\nBonferroni-Dunn vs control = {CONTROL} (rank diff > CD_BD={cd_bd:.2f} => sig):")
    rc = ranks[CONTROL]
    for m in ranks.index:
        if m == CONTROL: continue
        diff = ranks[m] - rc
        mark = "  *SIG*" if diff > cd_bd else ""
        log.append(f"   d_rank={diff:5.2f}  {CONTROL} vs {m:22s}{mark}")

    # 4. headline pairwise Wilcoxon (control vs each) + Holm
    others = [m for m in piv.columns if m != CONTROL]
    praw = []
    for m in others:
        try:
            _, pw = stats.wilcoxon(piv[CONTROL].values, piv[m].values)
        except ValueError:
            pw = 1.0
        praw.append(pw)
    order = np.argsort(praw)
    holm = np.empty(len(praw))
    for rank_i, idx in enumerate(order):
        holm[idx] = min(1.0, praw[idx] * (len(praw) - rank_i))
    holm = np.maximum.accumulate(holm[order])[np.argsort(order)]
    log.append(f"\nWilcoxon signed-rank, {CONTROL} vs each (Holm-adjusted):")
    for m, pr, ph in sorted(zip(others, praw, holm), key=lambda x: x[1]):
        log.append(f"   p_raw={pr:.3f}  p_holm={ph:.3f}  {'sig' if ph<ALPHA else 'ns '}  {m}")

    # CD diagram (Nemenyi)
    cd_diagram(ranks, cd_n,
               f"Nemenyi CD diagram - {tag}  (N={N}, k={k}, CD={cd_n:.2f})",
               os.path.join(OUT, f"cd_diagram_{split_mode}_{dataset}.png"))
    return piv


if __name__ == "__main__":
    df = pd.read_csv(SRC)
    log = ["Statistical comparison of classifier MCC (Demsar 2006 procedure)"]
    for ds in [d for d in ("full", "balanced") if d in set(df.dataset)]:     # the balanced arm is optional
        analyze(df, "cluster", ds, log)
    text = "\n".join(log)
    print(text)
    with open(os.path.join(OUT, "stats_report.txt"), "w") as fh:
        fh.write(text + "\n")
    print(f"\n[written] {OUT}/stats_report.txt  + nemenyi_*.csv + cd_diagram_*.png")
