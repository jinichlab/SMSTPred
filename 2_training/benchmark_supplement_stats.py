#!/usr/bin/env python3
"""Statistics quoted in the manuscript that are not produced by the plotting scripts.

Reads the bake-off results (sam_classifier_comparison_detailed.csv) and the per-iteration split sizes (iteration_splits.csv, parsed
from the bake-off log) and writes, into <bake-off output>/stats/:

  table_S4_wilcoxon_holm.csv   each other method versus the deployed MLP_ProstT5: mean MCC, mean paired difference, iterations won,
                               Wilcoxon signed-rank p, Holm-adjusted p over the fifteen comparisons          (manuscript Table S4)
  friedman_fig2A.txt           Friedman test and Bonferroni-Dunn comparisons for the nine methods of Fig 2A
  friedman_figS4.txt           the same for the eight methods of Fig S4
  eleven_iterations.txt        Fig 2A methods re-analysed on the eleven iterations that do not put the giant cluster in the test fold
  per_iteration_summary.csv    train / test size and MCC of the deployed model for each of the sixteen cluster-disjoint iterations
  iteration_summary.txt        the split-size and per-iteration numbers quoted in the Methods
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
BASE = HERE / os.environ.get("BAKEOFF_OUT", "sam_classifier_comparison_output")
OUT = BASE / "stats"
OUT.mkdir(exist_ok=True)
CONTROL = "MLP_ProstT5"
FIG2A = ["MLP_ProstT5", "MLP_SaProt", "MLP_ProtT5", "MLP_ESM-IF", "MLP_ProtT5+ESM-IF", "MLP_ProstT5+ESM-IF",
         "MMseqs2_1NN", "Foldseek_1NN", "Profile_HMM"]
FIGS4 = ["MLP_ProstT5", "SVM_Linear", "SVM_RBF", "SVM_Polynomial", "SVM_Sigmoidal", "Logistic Regression", "XGBoost", "Random Forest"]
GIANT_TRAIN = 2000          # iterations that put the ~5,700-protein cluster in the test fold leave fewer training proteins than this
ALPHA = 0.05

det = pd.read_csv(BASE / "sam_classifier_comparison_detailed.csv")
mcc = (det[(det.split_mode == "cluster") & (det.dataset == "full")]
       .pivot_table(index="iteration", columns="classifier", values="mcc"))
sizes = pd.read_csv(BASE / "iteration_splits.csv")
sizes = sizes[sizes.split_mode == "cluster"].set_index("iteration")


def holm(p):
    p = np.asarray(p, float); order = np.argsort(p); m = len(p)
    adj = np.empty(m); run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * p[i]); adj[i] = min(1.0, run)
    return adj


def friedman_block(df, methods, label, control=CONTROL):
    d = df[methods]; n, k = d.shape
    chi, p = stats.friedmanchisquare(*[d[m] for m in methods])
    ranks = d.apply(lambda r: pd.Series(stats.rankdata(-r.values), index=r.index), axis=1).mean().sort_values()
    cd = stats.norm.ppf(1 - ALPHA / (2 * (k - 1))) * np.sqrt(k * (k + 1) / (6 * n))
    lines = [f"{label}: N = {n} iterations, k = {k} methods",
             f"Friedman chi-square({k - 1}) = {chi:.1f}, p = {p:.2g}",
             f"Bonferroni-Dunn critical difference vs {control} = {cd:.2f} rank units", "", "mean rank (1 = best), mean MCC, rank difference vs control:"]
    sig = []
    for m, r in ranks.items():
        diff = r - ranks[control]
        flag = ""
        if m != control and abs(diff) > cd:
            flag = "  significantly different"; sig.append(m)
        lines.append(f"  {m:<22s} {r:5.2f}   MCC {d[m].mean():.3f}   {diff:+5.2f}{flag}")
    lines += ["", "significantly different from " + control + ": " + (", ".join(sig) if sig else "none")]
    return "\n".join(lines) + "\n", chi, p


# ---- Table S4: every other method versus the deployed MLP_ProstT5 (Wilcoxon signed-rank, Holm over the 15 comparisons)
others = [c for c in mcc.columns if c != CONTROL]
rows = []
for m in others:
    diff = mcc[m] - mcc[CONTROL]
    rows.append(dict(method=m, mean_mcc=mcc[m].mean(), mean_diff_vs_prostt5=diff.mean(), iterations_better=int((diff > 0).sum()),
                     iterations_worse=int((diff < 0).sum()), wilcoxon_p=stats.wilcoxon(mcc[m], mcc[CONTROL]).pvalue))
t4 = pd.DataFrame(rows)
t4["holm_p"] = holm(t4.wilcoxon_p)
t4["significant_holm_0.05"] = t4.holm_p < ALPHA
t4 = t4.sort_values("mean_mcc", ascending=False)
t4.insert(1, "n_iterations", len(mcc))
t4.to_csv(OUT / "table_S4_wilcoxon_holm.csv", index=False, float_format="%.6g")

# ---- Fig 2A and Fig S4
txt, chi, p = friedman_block(mcc, FIG2A, "Fig 2A")
(OUT / "friedman_fig2A.txt").write_text(txt)
txt, chi4, p4 = friedman_block(mcc, FIGS4, "Fig S4")
(OUT / "friedman_figS4.txt").write_text(txt)

# ---- eleven iterations with ordinary test sets
ordinary = sizes.index[sizes.n_train >= GIANT_TRAIN]
txt, chi11, p11 = friedman_block(mcc.loc[ordinary], FIG2A, "Fig 2A methods, eleven iterations with ordinary test sets")
w = stats.wilcoxon(mcc.loc[ordinary, "MLP_SaProt"], mcc.loc[ordinary, CONTROL])
txt += f"\nWilcoxon signed-rank, MLP_SaProt vs {CONTROL} on these iterations: p = {w.pvalue:.2f}\n"
(OUT / "eleven_iterations.txt").write_text(txt)

# ---- per-iteration table and the split-size numbers quoted in the Methods
per = sizes.join(mcc[[CONTROL]].rename(columns={CONTROL: "MLP_ProstT5_mcc"}))
per["giant_cluster_in_test"] = per.n_train < GIANT_TRAIN
per.to_csv(OUT / "per_iteration_summary.csv", float_format="%.6g")
g = per[per.giant_cluster_in_test]; o = per[~per.giant_cluster_in_test]
rand = pd.read_csv(BASE / "iteration_splits.csv").query("split_mode == 'random'")
lines = [f"cluster-disjoint iterations: {len(per)}",
         f"test-set size: {per.n_test.min()} to {per.n_test.max()}, median {per.n_test.median():g} (mean {per.n_test.mean():.0f}); random splits: constant {sorted(rand.n_test.unique())}",
         f"iterations with the giant cluster in the test fold: {', '.join(map(str, g.index))} ({len(g)} of {len(per)}); training proteins left: {g.n_train.min()} to {g.n_train.max()}",
         f"{CONTROL} MCC in those iterations: {g.MLP_ProstT5_mcc.min():.2f} to {g.MLP_ProstT5_mcc.max():.2f}; the other {len(o)}: {o.MLP_ProstT5_mcc.min():.2f} to {o.MLP_ProstT5_mcc.max():.2f}",
         f"{CONTROL} mean MCC: all {len(per)} = {per.MLP_ProstT5_mcc.mean():.3f}; the {len(o)} ordinary iterations = {o.MLP_ProstT5_mcc.mean():.3f}",
         f"Friedman on the ordinary iterations (Fig 2A methods): chi-square({len(FIG2A) - 1}) = {chi11:.1f}, p = {p11:.2g}"]
(OUT / "iteration_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"Fig 2A: chi2={chi:.1f} p={p:.2g} | Fig S4: chi2={chi4:.1f} p={p4:.2g}")
