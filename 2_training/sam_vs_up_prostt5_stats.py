#!/usr/bin/env python3
"""
Paired statistical comparison: SAM-trained vs UP-trained ProstT5 MLP, both
scored on the SAME held-out cluster-disjoint SAM test set, over N=16 repeated
90/10 splits (matches Panel A's design).

Because both models see identical test sets each iteration, this is a PAIRED
comparison over N=16 splits, per label + Overall (= mean across labels within
each iteration/model).

Tests (directional H1: SAM >= UP, focused training helps):
  - Paired t-test (naive; assumes independent iterations).
  - Wilcoxon signed-rank (nonparametric).
  - Nadeau-Bengio CORRECTED resampled t-test: inflates the variance by
    (1/n + n_test/n_train) to account for train-set overlap across repeated
    subsampling -> the statistically honest test here (this is the setting NB
    was originally derived for). rho = n_test/n_train = 0.1/0.9 = 1/9.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, os.environ.get("SAM_VS_UP_OUT", "sam_vs_up_prostt5_comparison_output"))
SRC = os.path.join(OUT, "sam_vs_up_prostt5_comparison_detailed.csv")
LABELS = ["Other", "Protein", "RNA", "DNA"]
N_ITER = 16
TRAIN_RATIO = 0.9
RHO = (1.0 - TRAIN_RATIO) / TRAIN_RATIO   # test/train ratio for repeated subsampling (Nadeau-Bengio)
METRIC = "mcc"


def corrected_t(diffs, rho=RHO):
    """Nadeau-Bengio corrected resampled paired t-test (one- and two-sided)."""
    n = len(diffs)
    m = np.mean(diffs)
    var = np.var(diffs, ddof=1)
    if var == 0:
        return m, np.inf, 0.0, 0.0
    se = np.sqrt((1.0 / n + rho) * var)
    t = m / se
    df = n - 1
    p_two = 2 * stats.t.sf(abs(t), df)
    p_one = stats.t.sf(t, df)      # H1: SAM > UP
    return m, t, p_one, p_two


def paired_table(d, metric):
    rows = []
    # per-label + Overall (mean across labels within each iteration/model)
    piv = d.pivot_table(index=["model", "iteration"], columns="label", values=metric)
    piv["Overall"] = piv[LABELS].mean(axis=1)
    for lab in LABELS + ["Overall"]:
        sam = piv.xs("SAM", level="model")[lab].sort_index().values
        up = piv.xs("UP", level="model")[lab].sort_index().values
        diffs = sam - up                       # >0 => SAM better
        # naive paired t (one-sided SAM>UP)
        t_naive, p_naive_two = stats.ttest_rel(sam, up)
        p_naive_one = p_naive_two / 2 if t_naive > 0 else 1 - p_naive_two / 2
        # Wilcoxon signed-rank (one-sided greater)
        try:
            w, p_wil = stats.wilcoxon(sam, up, alternative="greater")
        except ValueError:
            p_wil = np.nan
        # Nadeau-Bengio corrected
        m, t_c, p_c_one, p_c_two = corrected_t(diffs)
        rows.append(dict(
            label=lab, SAM=sam.mean(), UP=up.mean(), diff=m,
            t_naive=t_naive, p_naive_1sided=p_naive_one,
            p_wilcoxon_1sided=p_wil,
            t_corrected=t_c, p_corrected_1sided=p_c_one, p_corrected_2sided=p_c_two))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    d = pd.read_csv(SRC)
    log = [f"SAM vs UP ProstT5 MLP -- paired comparison over N={N_ITER} cluster-disjoint splits",
           f"Metric = {METRIC.upper()};  directional H1: SAM >= UP;  "
           f"Nadeau-Bengio rho = n_test/n_train = {RHO:.3f}", ""]

    tab = paired_table(d, METRIC)
    fmt = tab.copy()
    for c in fmt.columns:
        if c != "label":
            fmt[c] = fmt[c].map(lambda v: f"{v:.4f}" if pd.notna(v) else "nan")
    log.append(fmt.to_string(index=False))

    log.append("\nPer-iteration MCC (to show the variance driving the DNA result):")
    piv = d[d.label == "DNA"].pivot_table(index="iteration", columns="model", values="mcc")
    log.append("DNA per iteration:\n" + piv.round(3).to_string())

    log.append("\nReading:")
    log.append(" - p_corrected_1sided is the honest test (accounts for train-overlap across splits).")
    log.append(" - Nadeau-Bengio rho = n_test/n_train (repeated-subsampling setting).")
    text = "\n".join(log)
    print(text)
    tab.to_csv(os.path.join(OUT, "sam_vs_up_prostt5_stats.csv"), index=False)
    with open(os.path.join(OUT, "sam_vs_up_prostt5_stats.txt"), "w") as fh:
        fh.write(text + "\n")
    print(f"\n[written] {OUT}/sam_vs_up_prostt5_stats.csv + .txt")
