#!/usr/bin/env python3
"""Zoomed cluster-split MCC bar charts with mean labels and post-hoc significance.

Two method families, each drawn on the full and the class-balanced training set:

  Representations & tools (main Fig 2A) -- the six embedding MLPs vs the three
  off-the-shelf search-tool baselines (MMseqs2 sequence-homology, Profile HMM,
  Foldseek structure-homology). Asks: which INPUT REPRESENTATION wins?
    sam_classifier_cluster_mcc_zoom.png            full      (main Fig 2A)
    sam_classifier_cluster_mcc_zoom_balanced.png   balanced  (supplementary)

  MLP vs classical ML (supplementary) -- the deployed MLP_ProstT5 vs XGBoost /
  Random Forest / Logistic Regression / four SVMs, ALL trained on identical
  ProstT5 features. A controlled "same representation, different classifier head"
  contrast. Asks: given ProstT5, does the CLASSIFIER choice matter?
    sam_classifier_classical_ml_comparison.png            full      (supplementary)
    sam_classifier_classical_ml_comparison_balanced.png   balanced  (supplementary)

Per figure: cluster-disjoint MCC for each classifier sorted high->low, +/-1 SD
error bars, each bar labelled with its mean. Y-axis zoomed to amplify the top.

Statistics (Demsar 2006, from sam_classifier_comparison_detailed.csv), recomputed
WITHIN each figure's method set (so the critical difference matches k = methods
shown, not the full 16):
  - Friedman omnibus across the shown methods over the CV iterations.
  - Bonferroni-Dunn post-hoc vs control = MLP_ProstT5 (our deployment model):
    bars significantly WORSE than control (a=0.05) are starred/red; the leading
    run that is NOT significantly worse is bracketed as the "top tier".
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parent / "sam_classifier_comparison_output"
s = pd.read_csv(OUT / "sam_classifier_comparison_summary.csv")
d = pd.read_csv(OUT / "sam_classifier_comparison_detailed.csv")

CONTROL = "MLP_ProstT5"
ALPHA = 0.05
YMIN, YMAX = 0.50, 1.10
NITERS = int(d[d.split_mode == "cluster"].iteration.nunique())

# Method sets. CONTROL (MLP_ProstT5) is in both -- it anchors the classical-ML
# figure and is the Bonferroni-Dunn control everywhere.
EMBED_MLPS = {"MLP_ProstT5", "MLP_SaProt", "MLP_ProtT5", "MLP_ESM-IF",
              "MLP_ProtT5+ESM-IF", "MLP_ProstT5+ESM-IF"}
SEARCH_TOOLS = {"MMseqs2_1NN", "Profile_HMM", "Foldseek_1NN"}  # seq / profile / structure homology
CLASSICAL_ML = {"XGBoost", "Random Forest", "Logistic Regression",
                "SVM_Linear", "SVM_RBF", "SVM_Polynomial", "SVM_Sigmoidal"}

REPRO_TOOLS = EMBED_MLPS | SEARCH_TOOLS                 # Fig 2A method set
CLASSICAL_SET = CLASSICAL_ML | {CONTROL}               # supplementary method set


def posthoc(split_mode, dataset, methods):
    """Friedman p + set of methods significantly worse than CONTROL
    (Bonferroni-Dunn), from the paired iteration x method MCC matrix, restricted
    to `methods`."""
    sub = d[(d.split_mode == split_mode) & (d.dataset == dataset)
            & (d.classifier.isin(methods))]
    piv = sub.pivot_table(index="iteration", columns="classifier", values="mcc")
    N, k = piv.shape
    chi2, p = stats.friedmanchisquare(*[piv[c].values for c in piv.columns])
    ranks = piv.rank(axis=1, ascending=False).mean(axis=0)   # low = better
    z = stats.norm.ppf(1 - ALPHA / (2 * (k - 1)))
    cd = z * np.sqrt(k * (k + 1) / (6.0 * N))
    worse = {m for m in ranks.index
             if m != CONTROL and (ranks[m] - ranks[CONTROL]) > cd}
    return dict(chi2=chi2, p=p, df=k - 1, N=N, cd=cd, worse=worse)


def color(clf):
    if clf in EMBED_MLPS:
        return "#2c7fb8"          # embedding MLPs
    if clf in SEARCH_TOOLS:
        return "#e6550d"          # homology / structure-search baselines
    return "#969696"              # classical ML


def draw_panel(ax, ds, methods, bracket_label=f"top tier — n.s. vs {CONTROL}"):
    sub = (s[(s.split_mode == "cluster") & (s.dataset == ds)
             & (s.classifier.isin(methods))]
           .sort_values("mcc_mean", ascending=False).reset_index(drop=True))
    st = posthoc("cluster", ds, methods)
    x = range(len(sub))
    ax.bar(x, sub.mcc_mean, yerr=sub.mcc_std, capsize=3,
           color=[color(c) for c in sub.classifier], edgecolor="black",
           linewidth=0.6, error_kw=dict(ecolor="#444444", lw=1))
    for i, (mm, sdd) in enumerate(zip(sub.mcc_mean, sub.mcc_std)):
        ax.text(i, min(mm + sdd + 0.006, 0.998), f"{mm:.3f}",
                ha="center", va="bottom", fontsize=9, fontweight="bold")
    # top-tier band + bracket: leading run NOT significantly worse than control
    last = 0
    for i, c in enumerate(sub.classifier):
        if c in st["worse"]:
            break
        last = i
    ax.axvspan(-0.5, last + 0.5, color="#08519c", alpha=0.06, zorder=0)
    yb = 1.02
    ax.plot([-0.4, last + 0.4], [yb, yb], color="#08519c", lw=1.4, clip_on=False)
    for xe in (-0.4, last + 0.4):
        ax.plot([xe, xe], [yb - 0.012, yb], color="#08519c", lw=1.4, clip_on=False)
    ax.text(last / 2.0, yb + 0.008, bracket_label,
            ha="center", va="bottom", fontsize=9.5, fontweight="bold",
            color="#08519c", clip_on=False)
    # x tick labels: control bold-blue; significantly worse = red + trailing *
    ax.set_xticks(list(x))
    ticklabels = [f"{c} *" if c in st["worse"] else c for c in sub.classifier]
    ax.set_xticklabels(ticklabels, rotation=45, ha="right", fontsize=9)
    for lab, clf in zip(ax.get_xticklabels(), sub.classifier):
        if clf == CONTROL:
            lab.set_color("#08519c"); lab.set_fontweight("bold")
        elif clf in st["worse"]:
            lab.set_color("crimson")
    ax.set_ylim(YMIN, YMAX)
    ax.set_yticks(np.arange(0.5, 1.001, 0.1))
    ax.set_ylabel(f"MCC (mean ± SD over {NITERS} iters)", fontsize=11)
    ax.grid(True, axis="y", alpha=0.3)
    ax.axhline(YMIN, color="black", lw=0.8)
    ax.text(0.985, 0.90,
            f"Friedman $\\chi^2$({st['df']})={st['chi2']:.1f}, p={st['p']:.1e}\n"
            f"(N={st['N']} iters, methods differ overall)",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.35", fc="#f2f2f2", ec="#999", lw=0.6))
    return st


# Legend used by plot_fig2_main for the composite (Fig 2A = representations & tools).
LEGEND = [mpatches.Patch(color="#2c7fb8", label="Embedding MLP"),
          mpatches.Patch(color="#e6550d", label="Homology / structure-search baseline")]

LEGEND_CLASSICAL = [mpatches.Patch(color="#2c7fb8", label="Embedding MLP (deployed)"),
                    mpatches.Patch(color="#969696", label="Classical ML (on ProstT5)")]


def _saprot_prostt5_note(fig, ds):
    """Bottom caption line: SaProt vs ProstT5 tie (only when both are shown)."""
    pv = (d[(d.split_mode == "cluster") & (d.dataset == ds)]
          .pivot_table(index="iteration", columns="classifier", values="mcc"))
    dd = pv["MLP_SaProt"] - pv["MLP_ProstT5"]
    wp = stats.wilcoxon(pv["MLP_SaProt"], pv["MLP_ProstT5"]).pvalue
    fig.text(0.5, 0.060,
             "MLP_SaProt and MLP_ProstT5 are the top-ranked structure-aware MLPs and are statistically tied",
             ha="center", va="bottom", fontsize=8.5, color="#333")
    fig.text(0.5, 0.035,
             f"(Δ={dd.mean():+.3f} MCC, SaProt wins {int((dd > 0).sum())}/{len(dd)} iters, Wilcoxon p={wp:.2f}); "
             f"ProstT5 is deployed as it needs only sequence input.",
             ha="center", va="bottom", fontsize=8.5, color="#333")


def make_fig(ds, methods, fname, legend, bracket_label, saprot_note):
    fig, ax = plt.subplots(figsize=(12.5, 7.0))
    draw_panel(ax, ds, methods, bracket_label=bracket_label)
    fig.legend(handles=legend, fontsize=9.5, loc="upper center", ncol=len(legend),
               bbox_to_anchor=(0.5, 0.985), frameon=True)
    if saprot_note:
        _saprot_prostt5_note(fig, ds)
    fig.tight_layout(rect=[0, 0.10 if saprot_note else 0.03, 1, 0.95])
    p = OUT / fname
    fig.savefig(p, dpi=150)
    plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    # --- Fig 2A + balanced supplementary: representations & tools ---
    make_fig("full", REPRO_TOOLS, "sam_classifier_cluster_mcc_zoom.png",
             LEGEND, f"top tier — n.s. vs {CONTROL}", saprot_note=True)
    make_fig("balanced", REPRO_TOOLS, "sam_classifier_cluster_mcc_zoom_balanced.png",
             LEGEND, f"top tier — n.s. vs {CONTROL}", saprot_note=True)

    # --- supplementary: MLP vs classical ML, all on ProstT5 ---
    make_fig("full", CLASSICAL_SET, "sam_classifier_classical_ml_comparison.png",
             LEGEND_CLASSICAL, f"top tier — n.s. vs {CONTROL}", saprot_note=False)
    make_fig("balanced", CLASSICAL_SET, "sam_classifier_classical_ml_comparison_balanced.png",
             LEGEND_CLASSICAL, f"top tier — n.s. vs {CONTROL}", saprot_note=False)

    for name, methods in [("representations+tools", REPRO_TOOLS), ("classical-ML", CLASSICAL_SET)]:
        for ds in ("full", "balanced"):
            st = posthoc("cluster", ds, methods)
            print(f"  {name} cluster/{ds}: Friedman p={st['p']:.2e}; "
                  f"worse-than-{CONTROL} = {sorted(st['worse'])}")
