#!/usr/bin/env python3
"""Fig 5: SAM-methyltransferase counts against proteome size in eukaryotes, one fit per clade (plants, fungi, metazoa, other eukaryotes).

  python scaling_by_clade.py --dir <release>/fig5        (after build_proteome_vs_sam_table.py; run_fig5.sh calls it)

Why by clade: the single "other eukaryotes" group of scaling_loglog.py panels C and D mixes clades that sit at very different levels (metabolite-substrate MTases at
20,000 proteins: about 100 in plants and fungi, about 12 in metazoa), so one power law through it fits poorly (reduced chi-square 29 and 54). Here every clade is fitted
separately with the same binning and weighted non-linear least squares as scaling_loglog.py (0.15-decade size bins, at least MIN_BIN proteomes per bin, n = lambda g^alpha),
with 95% intervals from clade-level bootstrap resampling (lineage prefix of depth 8).

Outputs (in --dir): fig5_by_clade.png/.pdf, fig5_by_clade.csv (per clade and measure: n, alpha with interval, usable bins, reduced chi-square, fitted height at 20,000 proteins with
interval) and fig5_clade_model_comparison.csv (log-log least squares on the proteomes: one line, clade offsets with a common slope, clade-specific lines; R2, RMSE, AIC, BIC).
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import scaling_loglog as sl     # noqa: E402

REF_SIZE = 2e4                  # proteins; the size at which the fitted heights are compared
N_BOOT = 1000
CLADES = ["Metazoa", "Fungi", "Other eukaryotes", "Plants"]
COL = {"Plants": "#618d43", "Fungi": "#8c5a2b", "Metazoa": "#3a6fb0", "Other eukaryotes": "#9b6bb5"}
MEASURES = [("n_sam", "Total SAM MTases"), ("n_Metabolite", "Metabolite-substrate SAM MTases")]


def clade_of(lineage):
    for key, name in (("Viridiplantae", "Plants"), ("Metazoa", "Metazoa"), ("Fungi", "Fungi")):
        if key in lineage:
            return name
    return "Other eukaryotes"


def compare_models(d, col):
    """Log-log least squares on the proteomes: one line, clade offsets with a common slope, clade-specific lines."""
    x, y, cl = np.log10(d.proteome_size.values.astype(float)), np.log10(d[col].values.astype(float)), d.clade.values
    n = len(d)
    levels = [c for c in CLADES if c in set(cl)]
    D = np.column_stack([(cl == c).astype(float) for c in levels[1:]])
    designs = {"one line": np.c_[np.ones(n), x], "clade offsets, common slope": np.c_[np.ones(n), x, D], "clade-specific lines": np.c_[np.ones(n), x, D, D * x[:, None]]}
    rows = []
    for nm, X in designs.items():
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        rss = float(((y - X @ b) ** 2).sum()); k = X.shape[1]
        rows.append(dict(measure=col, model=nm, n=n, parameters=k, R2=1 - rss / float(((y - y.mean()) ** 2).sum()), RMSE_log10=float(np.sqrt(rss / n)),
                         AIC=n * np.log(rss / n) + 2 * k, BIC=n * np.log(rss / n) + k * np.log(n)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="the release's fig5 folder (holds proteome_vs_sam_bylabel.csv)")
    a = ap.parse_args()
    out = Path(a.dir)
    sl.IO = out
    t = sl.build_table()
    e = t[t.domain == "Eukaryota"].copy()
    e["clade"] = e.lineage.map(clade_of)
    ed = sl.edges()
    rng = np.random.default_rng(sl.SEED)
    rows, comp = [], []
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))
    for ax, (col, ylab), letter in zip(axes, MEASURES, "AB"):
        d_all = e[e[col] > 0]
        comp += compare_models(d_all, col)
        for c in CLADES:
            d = d_all[d_all.clade == c]
            g, y = d.proteome_size.values.astype(float), d[col].values.astype(float)
            ax.scatter(g, y, s=14, color=COL[c], alpha=0.35, edgecolor="none")
            n, mg, my, se, u = sl.bins(g, y, ed, sl.MIN_BIN)
            ax.errorbar(mg[u], my[u], yerr=se[u], fmt="s", ms=6, mfc="white", mec=COL[c], color=COL[c], lw=1.2)
            al, la, chi, nb = sl.wnls(g, y, ed)
            height = float(np.exp(la) * (REF_SIZE / sl.G0) ** al) if np.isfinite(al) else np.nan
            groups = sl.groups_of(d[f"clade{sl.CLUSTER_DEPTH['Eukaryota']}"].values)
            bs = []
            for _ in range(N_BOOT):
                ix = sl.resample(groups, rng)
                r = sl.wnls(g[ix], y[ix], ed)
                if np.isfinite(r[0]):
                    bs.append((r[0], np.exp(r[1]) * (REF_SIZE / sl.G0) ** r[0]))
            bs = np.array(bs) if len(bs) > 100 else np.full((2, 2), np.nan)
            alo, ahi = np.nanpercentile(bs[:, 0], [2.5, 97.5]); hlo, hhi = np.nanpercentile(bs[:, 1], [2.5, 97.5])
            if np.isfinite(al):
                xs = np.geomspace(mg[u].min(), mg[u].max(), 50)
                ax.plot(xs, np.exp(la) * (xs / sl.G0) ** al, color=COL[c], lw=2.2, label=f"{c} (n={len(d)}): α = {al:.2f} [{alo:.2f}, {ahi:.2f}]")
            else:
                ax.plot([], [], color=COL[c], lw=2.2, label=f"{c} (n={len(d)}): not estimable")
            rows.append(dict(measure=col, clade=c, n=len(d), alpha=al, alpha_lo=alo, alpha_hi=ahi, usable_bins=nb, reduced_chi2=chi, height_at_20k=height, height_lo=hlo, height_hi=hhi))
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Proteome size (number of proteins)", fontsize=11); ax.set_ylabel(ylab, fontsize=11)
        ax.grid(True, which="major", ls=":", alpha=0.35); ax.axvline(REF_SIZE, color="k", lw=0.6, ls=":")
        ax.legend(loc="upper left", fontsize=9, frameon=True)
        ax.text(0.0, 1.02, letter, transform=ax.transAxes, fontsize=20, fontweight="bold", ha="left", va="bottom")
    fig.tight_layout()
    fig.savefig(out / "fig5_by_clade.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / "fig5_by_clade.pdf", bbox_inches="tight")
    pd.DataFrame(rows).to_csv(out / "fig5_by_clade.csv", index=False, float_format="%.4g")
    pd.DataFrame(comp).to_csv(out / "fig5_clade_model_comparison.csv", index=False, float_format="%.5g")
    print("wrote fig5_by_clade.png/.pdf, fig5_by_clade.csv, fig5_clade_model_comparison.csv")
    print(pd.DataFrame(rows).round(2).to_string(index=False))
    print(pd.DataFrame(comp).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
