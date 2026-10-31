#!/usr/bin/env python3
"""Zoomed cluster-split MCC bar charts with mean labels and post-hoc significance.

Two method families, drawn on the full training set:

  Representations & tools (main Fig 2A) -- the six embedding MLPs vs the three
  off-the-shelf search-tool baselines (MMseqs2 sequence-homology, Profile HMM,
  Foldseek structure-homology). Asks: which INPUT REPRESENTATION wins?
    sam_classifier_cluster_mcc_zoom.png

  MLP vs classical ML (Fig S4) -- the deployed MLP_ProstT5 vs XGBoost /
  Random Forest / Logistic Regression / four SVMs, ALL trained on identical
  ProstT5 features. A controlled "same representation, different classifier head"
  contrast. Asks: given ProstT5, does the CLASSIFIER choice matter?
    sam_classifier_classical_ml_comparison.png

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
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parent / os.environ.get("BAKEOFF_OUT", "sam_classifier_comparison_output")   # BAKEOFF_OUT: optional re-run folder
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


# Style knobs. The defaults are the look of the supplementary figures (Fig S4 and the standalone Fig 2A); plot_fig2_main passes its own
# (softer colours, larger type) so the composite Figure 2 can differ without changing the others.
STYLE = dict(embed="#4578aa", search="#ed8c44", other="#969696", accent="#2e5a86",
             fs_tick=9, fs_ylabel=11, fs_value=9, fs_bracket=9.5, fs_stats=9, ytick_pad=3.5,
             value_bbox=False, stats_short=False, stats_xy=(0.985, 0.90))


def color(clf, style=STYLE):
    if clf in EMBED_MLPS:
        return style["embed"]         # embedding MLPs
    if clf in SEARCH_TOOLS:
        return style["search"]        # homology / structure-search baselines
    return style["other"]             # classical ML


def draw_panel(ax, ds, methods, bracket_label=f"top tier — n.s. vs {CONTROL}", style=STYLE):
    sub = (s[(s.split_mode == "cluster") & (s.dataset == ds)
             & (s.classifier.isin(methods))]
           .sort_values("mcc_mean", ascending=False).reset_index(drop=True))
    st = posthoc("cluster", ds, methods)
    x = range(len(sub))
    ax.bar(x, sub.mcc_mean, yerr=sub.mcc_std, capsize=3,
           color=[color(c, style) for c in sub.classifier], edgecolor="black",
           linewidth=0.6, error_kw=dict(ecolor="#444444", lw=1))
    for i, (mm, sdd) in enumerate(zip(sub.mcc_mean, sub.mcc_std)):
        ax.text(i, min(mm + sdd + 0.006, 0.998), f"{mm:.3f}",
                ha="center", va="bottom", fontsize=style["fs_value"], fontweight="bold",
                bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.9) if style["value_bbox"] else None)
    # top-tier band + bracket: leading run NOT significantly worse than control
    last = 0
    for i, c in enumerate(sub.classifier):
        if c in st["worse"]:
            break
        last = i
    ax.axvspan(-0.5, last + 0.5, color=style["accent"], alpha=0.06, zorder=0)
    yb = 1.02
    ax.plot([-0.4, last + 0.4], [yb, yb], color=style["accent"], lw=1.4, clip_on=False)
    for xe in (-0.4, last + 0.4):
        ax.plot([xe, xe], [yb - 0.012, yb], color=style["accent"], lw=1.4, clip_on=False)
    ax.text(last / 2.0, yb + 0.008, bracket_label,
            ha="center", va="bottom", fontsize=style["fs_bracket"], fontweight="bold",
            color=style["accent"], clip_on=False)
    # x tick labels: control bold-blue; significantly worse = red + trailing *
    ax.set_xticks(list(x))
    ticklabels = [f"{c} *" if c in st["worse"] else c for c in sub.classifier]
    ax.set_xticklabels(ticklabels, rotation=45, ha="right", fontsize=style["fs_tick"])
    for lab, clf in zip(ax.get_xticklabels(), sub.classifier):
        if clf == CONTROL:
            lab.set_color(style["accent"]); lab.set_fontweight("bold")
        elif clf in st["worse"]:
            lab.set_color("crimson")
    ax.set_ylim(YMIN, YMAX)
    ax.set_yticks(np.arange(0.5, 1.001, 0.1))
    ax.set_ylabel(f"MCC (mean ± SD over {NITERS} iters)", fontsize=style["fs_ylabel"])
    ax.tick_params(axis="y", labelsize=style["fs_tick"], pad=style["ytick_pad"])
    ax.grid(True, axis="y", alpha=0.3)
    ax.axhline(YMIN, color="black", lw=0.8)
    stats_txt = (f"Friedman $\\chi^2$({st['df']}) = {st['chi2']:.1f}, p = {st['p']:.0e}" if style["stats_short"] else
                 f"Friedman $\\chi^2$({st['df']})={st['chi2']:.1f}, p={st['p']:.1e}\n"
                 f"(N={st['N']} iters, methods differ overall)")
    ax.text(*style["stats_xy"], stats_txt,
            transform=ax.transAxes, ha="right", va="top", fontsize=style["fs_stats"],
            bbox=dict(boxstyle="round,pad=0.35", fc="#f2f2f2", ec="#999", lw=0.6))
    return st


# Legend used by plot_fig2_main for the composite (Fig 2A = representations & tools).
def legend_handles(style=STYLE):
    return [mpatches.Patch(facecolor=style["embed"], edgecolor="black", linewidth=0.6, label="Embedding MLP"),
            mpatches.Patch(facecolor=style["search"], edgecolor="black", linewidth=0.6, label="Homology / structure-search baseline")]


LEGEND = [mpatches.Patch(color="#4578aa", label="Embedding MLP"),
          mpatches.Patch(color="#ed8c44", label="Homology / structure-search baseline")]

LEGEND_CLASSICAL = [mpatches.Patch(color="#4578aa", label="Embedding MLP (deployed)"),
                    mpatches.Patch(color="#969696", label="Classical ML (on ProstT5)")]


def _saprot_prostt5_note(fig, ds):
    """Bottom caption line: deployed ProstT5 (sequence only) vs SaProt (needs structure)."""
    pv = (d[(d.split_mode == "cluster") & (d.dataset == ds)]
          .pivot_table(index="iteration", columns="classifier", values="mcc"))
    dd = pv["MLP_SaProt"] - pv["MLP_ProstT5"]
    wp = stats.wilcoxon(pv["MLP_SaProt"], pv["MLP_ProstT5"]).pvalue
    verdict = "statistically tied" if wp >= ALPHA else "significantly different"
    fig.text(0.5, 0.060,
             f"MLP_ProstT5 (deployed; needs only sequence) and MLP_SaProt (needs a structure) are {verdict}",
             ha="center", va="bottom", fontsize=8.5, color="#333")
    fig.text(0.5, 0.035,
             f"(SaProt minus ProstT5: Δ={dd.mean():+.3f} MCC, SaProt wins {int((dd > 0).sum())}/{len(dd)} iters, Wilcoxon p={wp:.2f}).",
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
    # --- Fig 2A: representations & tools ---
    make_fig("full", REPRO_TOOLS, "sam_classifier_cluster_mcc_zoom.png",
             LEGEND, f"top tier — n.s. vs {CONTROL}", saprot_note=True)

    # --- Fig S4: MLP vs classical ML, all on ProstT5 ---
    make_fig("full", CLASSICAL_SET, "sam_classifier_classical_ml_comparison.png",
             LEGEND_CLASSICAL, f"top tier — n.s. vs {CONTROL}", saprot_note=False)

    for name, methods in [("representations+tools", REPRO_TOOLS), ("classical-ML", CLASSICAL_SET)]:
        st = posthoc("cluster", "full", methods)
        print(f"  {name} cluster/full: Friedman chi2({st['df']})={st['chi2']:.1f}, p={st['p']:.2e}; "
              f"worse-than-{CONTROL} = {sorted(st['worse'])}")
