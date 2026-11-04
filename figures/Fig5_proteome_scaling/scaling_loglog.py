#!/usr/bin/env python3
"""Fig 5 (log-log): SAM-MTase count vs proteome size, lineage vs rest of its domain, as power laws n = lambda * g^alpha.

Adapted for a refresh release from the scaling-laws-fig5 branch (analysis by A. Jinich). Panels: A Bacteria total, B Bacteria metabolite
(Actinomycetota vs other bacteria), C Eukaryota total, D Eukaryota metabolite (Viridiplantae vs other eukaryotes).

  python scaling_loglog.py [--dir <release>/fig5]       (default: <newest release>/fig5; needs proteome_vs_sam_bylabel.csv from build_proteome_vs_sam_table.py)

Input : <release>/fig5/proteome_vs_sam_bylabel.csv (per-proteome size, IPR029063 count n_sam, ensemble calls n_<label>) + lineage from data/selected_1259_organisms.csv
Output: <dir>/fig5_scaling_loglog.png/.pdf, fig5_scaling_numbers.csv (every exponent, CI, test), fig5_scaling_bins.csv, fig5_scaling_table.csv

Fixed analysis settings (one choice, see FIG5_SCALING_NOTES.md): proteomes log-binned by size (0.15 decades from 10^3), per bin the mean count and its
standard error, bins with fewer than MIN_BIN = 5 proteomes dropped (10 as in the 2017 NAR protocol leaves Viridiplantae one bin; the strict version is
reported as analysis 'strict_protocol'), weighted non-linear least squares of n = exp(a) (g/1e4)^alpha. 95% CIs and the lineage-minus-rest difference come
from 2,000 bootstrap replicates that resample whole clades (lineage prefix of depth 7 in bacteria = family level, 8 in eukaryotes); lineage and rest are
resampled independently. The 'fold at the median lineage size' is the level difference at equal proteome size.
"""
import argparse
import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
IO = HERE
DATA = Path(os.environ.get("SMSTPRED_DATA", REPO / "data"))
META = next((p for p in (DATA / "selected_1259_organisms.csv", REPO / "data" / "selected_1259_organisms.csv") if p.exists()), None)
LABELS = ["DNA", "RNA", "Protein", "Metabolite"]

BIN_W, BIN_START, MIN_BIN, STRICT_MIN = 0.15, 3.0, 5, 10
G0 = 1e4
B = 2000
SEED = 20261102
CLUSTER_DEPTH = {"Bacteria": 7, "Eukaryota": 8}

BASE_COLOR, ACTINO_COLOR, PLANT_COLOR, FIT_COLOR = "#707070", "#ed8c44", "#618d43", "#222222"
PANELS = [  # letter, domain, count column, lineage column, lineage name, rest name, y label, colour
    ("A", "Bacteria", "n_sam", "is_actino", "Actinomycetota", "Other bacteria", "Total SAM MTases", ACTINO_COLOR),
    ("B", "Bacteria", "n_Metabolite", "is_actino", "Actinomycetota", "Other bacteria", "Metabolite-substrate SAM MTases", ACTINO_COLOR),
    ("C", "Eukaryota", "n_sam", "is_plant", "Viridiplantae", "Other eukaryotes", "Total SAM MTases", PLANT_COLOR),
    ("D", "Eukaryota", "n_Metabolite", "is_plant", "Viridiplantae", "Other eukaryotes", "Metabolite-substrate SAM MTases", PLANT_COLOR),
]


def canon(n):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9_\-]+", "", re.sub(r"\s+", "_", str(n).strip()))).strip("_")


# ----------------------------------------------------------------------------------------------------------------------
# 1. table
# ----------------------------------------------------------------------------------------------------------------------
def build_table():
    t = pd.read_csv(IO / "proteome_vs_sam_bylabel.csv")
    t = t[t.domain.isin(["Bacteria", "Eukaryota"])].reset_index(drop=True)      # archaea: too few proteomes to fit
    meta = pd.read_csv(META)
    meta["key"] = meta.Organism.map(canon)
    t["key"] = t.organism.map(canon)
    t = t.merge(meta[["key", "Taxonomic lineage", "genus"]].drop_duplicates("key"), on="key", how="left").drop(columns="key")
    t["lineage"] = t["Taxonomic lineage"].fillna("")
    t["genus"] = t.genus.fillna(t.organism.str.split("_").str[0])
    t["is_actino"] = t.lineage.str.contains("Actinomycetota")
    t["is_plant"] = t.lineage.str.contains("Viridiplantae")
    for d in (6, 7, 8):
        t[f"clade{d}"] = t.lineage.map(lambda s, d=d: ", ".join(s.split(", ")[:d]) if s else "")
        t.loc[t[f"clade{d}"] == "", f"clade{d}"] = t.genus          # the organism without a lineage string is its own clade
    print(f"table: {len(t)} proteomes ({t.domain.value_counts().to_dict()}), {int(t.n_sam.sum())} SAM-MTases", flush=True)
    return t.drop(columns=["Taxonomic lineage"])


# ----------------------------------------------------------------------------------------------------------------------
# 2. fit
# ----------------------------------------------------------------------------------------------------------------------
def edges(w=BIN_W):
    return 10 ** np.arange(BIN_START, 5.4 + w, w)


def bins(g, y, ed, minbin):
    k = np.digitize(g, ed) - 1
    nb = len(ed) - 1
    ok = (k >= 0) & (k < nb)
    k, g, y = k[ok], g[ok], y[ok]
    n = np.bincount(k, minlength=nb).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        mg, my = np.bincount(k, g, nb) / n, np.bincount(k, y, nb) / n
        se = np.sqrt(np.maximum((np.bincount(k, y * y, nb) - n * my ** 2) / (n - 1), 0) / n)
    used = (n >= minbin) & (se > 0) & (my > 0)
    return n, mg, my, se, used


def wnls(g, y, ed, minbin=MIN_BIN):
    """-> (alpha, a, chi2_red, nbins); NaNs when fewer than 3 usable bins."""
    n, mg, my, se, u = bins(g, y, ed, minbin)
    if u.sum() < 3:
        return np.nan, np.nan, np.nan, int(u.sum())
    x, m, s = mg[u], my[u], se[u]
    f = lambda xx, a, al: np.exp(a) * (xx / G0) ** al
    p0 = np.polyfit(np.log(x / G0), np.log(m), 1, w=m / s)
    try:
        p, _ = curve_fit(f, x, m, p0=[p0[1], p0[0]], sigma=s, absolute_sigma=True, maxfev=5000)
    except Exception:
        return np.nan, np.nan, np.nan, int(u.sum())
    return p[1], p[0], float(np.sum(((m - f(x, *p)) / s) ** 2) / max(u.sum() - 2, 1)), int(u.sum())


def pct(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return tuple(np.percentile(v, [2.5, 97.5])) if len(v) > 100 else (np.nan, np.nan)


def groups_of(labels):
    return [np.where(labels == k)[0] for k in pd.unique(labels)]


def resample(groups, rng):
    return np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])


def pval(v):
    v = v[np.isfinite(v)]
    return float(min(1.0, 2 * min((v <= 0).mean(), (v >= 0).mean()))) if len(v) else np.nan


def analyse(panel, t, rng, out, binrows):
    letter, dm, ycol, lcol, lname, rname, _, _ = panel
    s = t[t.domain == dm].reset_index(drop=True)
    L = s[lcol].values
    g, y = s.proteome_size.values.astype(float), s[ycol].values.astype(float)
    gl, yl, gr, yr = g[L], y[L], g[~L], y[~L]
    ed = edges()
    gref = float(np.median(gl))                                   # reference size for the offset: median size of the lineage
    fold = lambda al_l, a_l, al_r, a_r: np.exp(a_l + al_l * np.log(gref / G0) - a_r - al_r * np.log(gref / G0))
    row = lambda analysis, group, q, est, lo=np.nan, hi=np.nan, n=np.nan, note="": out.append(dict(
        panel=letter, domain=dm, measure=ycol, analysis=analysis, group=group, quantity=q, estimate=est, ci_lo=lo, ci_hi=hi, n=n, note=note))

    def pair_boot(clade_col, minbin=MIN_BIN, rest_mask=None):
        gg_l = groups_of(s.loc[L, clade_col].values)
        rm = ~L if rest_mask is None else (~L) & rest_mask
        gg_r = groups_of(s.loc[rm, clade_col].values)
        gr_, yr_ = g[rm], y[rm]
        res = []
        for _ in range(B):
            il, ir = resample(gg_l, rng), resample(gg_r, rng)
            al_l, a_l, _, _ = wnls(gl[il], yl[il], ed, minbin)
            al_r, a_r, _, _ = wnls(gr_[ir], yr_[ir], ed, minbin)
            res.append((al_l, al_r, fold(al_l, a_l, al_r, a_r)))
        return np.array(res)

    # ---- primary ---------------------------------------------------------------------------------------------------
    fits = {}
    for nm, gg_, yy_ in ((lname, gl, yl), (rname, gr, yr)):
        al, a, chi, nb = wnls(gg_, yy_, ed)
        fits[nm] = (al, a)
        n, mg, my, se, u = bins(gg_, yy_, ed, MIN_BIN)
        for i in np.where(n > 0)[0]:
            binrows.append(dict(panel=letter, domain=dm, measure=ycol, group=nm, bin_lo=ed[i], bin_hi=ed[i + 1], n_proteomes=int(n[i]),
                                mean_size=mg[i], mean_count=my[i], se_count=se[i], used_in_fit=bool(u[i])))
        row("primary", nm, "alpha", al, n=len(gg_), note=f"{nb} bins; lambda={np.exp(a) * G0 ** -al:.4g}")
        row("primary", nm, "reduced_chi2", chi, n=len(gg_), note=f"{nb} bins")
    bt = pair_boot(f"clade{CLUSTER_DEPTH[dm]}")
    for j, nm in enumerate((lname, rname)):
        lo, hi = pct(bt[:, j])
        out[[i for i, r in enumerate(out) if r["panel"] == letter and r["analysis"] == "primary" and r["group"] == nm and r["quantity"] == "alpha"][0]].update(ci_lo=lo, ci_hi=hi)
    (al_l, a_l), (al_r, a_r) = fits[lname], fits[rname]
    dlo, dhi = pct(bt[:, 0] - bt[:, 1])
    flo, fhi = pct(bt[:, 2])
    ci_l, ci_r = pct(bt[:, 0]), pct(bt[:, 1])
    row("primary", "lineage - rest", "delta_alpha", al_l - al_r, dlo, dhi, note=f"p_boot={pval(bt[:, 0] - bt[:, 1]):.3f}; individual CIs overlap: {not (ci_l[0] > ci_r[1] or ci_r[0] > ci_l[1])}")
    row("primary", "lineage / rest", "fold_at_median_lineage_size", float(fold(al_l, a_l, al_r, a_r)), flo, fhi, note=f"size={gref:.0f} proteins; p_boot(fold<=1)={pval(np.log(bt[:, 2])):.3f}")
    row("primary", "lineage / rest", "bootstrap_failures", float(np.isnan(bt[:, 0] - bt[:, 1]).mean()), note="share of replicates with <3 usable bins in either group")

    # ---- sensitivities ----------------------------------------------------------------------------------------------
    sl = [wnls(gl, yl, ed, STRICT_MIN), wnls(gr, yr, ed, STRICT_MIN)]
    for nm, f_ in zip((lname, rname), sl):
        row("strict_protocol", nm, "alpha", f_[0], n=f_[3], note=f"min {STRICT_MIN} proteomes per bin: {f_[3]} usable bins" + ("" if f_[3] >= 3 else " -> not estimable"))
    bt6 = pair_boot("clade6")
    lo, hi = pct(bt6[:, 0] - bt6[:, 1])
    row("coarser_clades_depth6", "lineage - rest", "delta_alpha", al_l - al_r, lo, hi, note=f"p_boot={pval(bt6[:, 0] - bt6[:, 1]):.3f}")
    row("coarser_clades_depth6", "lineage / rest", "fold_at_median_lineage_size", float(fold(al_l, a_l, al_r, a_r)), *pct(bt6[:, 2]))
    mm = (g >= gl.min()) & (g <= gl.max())
    al_rm, a_rm, _, nb_rm = wnls(g[(~L) & mm], y[(~L) & mm], ed)
    btm = pair_boot(f"clade{CLUSTER_DEPTH[dm]}", rest_mask=mm)
    lo, hi = pct(btm[:, 0] - btm[:, 1])
    row("rest_restricted_to_lineage_size_range", rname, "alpha", al_rm, *pct(btm[:, 1]), n=int(((~L) & mm).sum()), note=f"{nb_rm} bins")
    row("rest_restricted_to_lineage_size_range", "lineage - rest", "delta_alpha", al_l - al_rm, lo, hi, note=f"p_boot={pval(btm[:, 0] - btm[:, 1]):.3f}")
    # unbinned cross-check: Poisson GLM with lineage x log-size interaction, clade-clustered SEs; size centred at gref so 'lineage'
    # is the offset at gref.
    X = pd.DataFrame({"c": 1.0, "ls": np.log(g / gref), "L": L.astype(float)})
    X["Lls"] = X.L * X.ls
    r = sm.GLM(y, X, family=sm.families.Poisson()).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(s[f"clade{CLUSTER_DEPTH[dm]}"])[0]})
    ci = r.conf_int()
    row("poisson_glm", rname, "alpha", r.params["ls"], ci.loc["ls", 0], ci.loc["ls", 1], n=len(s))
    row("poisson_glm", lname, "alpha", r.params["ls"] + r.params["Lls"], note="rest + interaction")
    row("poisson_glm", "lineage - rest", "delta_alpha", r.params["Lls"], ci.loc["Lls", 0], ci.loc["Lls", 1], note=f"p={r.pvalues['Lls']:.3g}")
    row("poisson_glm", "lineage / rest", "fold_at_median_lineage_size", float(np.exp(r.params["L"])), float(np.exp(ci.loc["L", 0])), float(np.exp(ci.loc["L", 1])), note=f"p={r.pvalues['L']:.3g}")
    # the previous Fig 5 statistic: ratio of OLS slopes of untransformed counts
    ols = lambda m: np.polyfit(g[m], y[m], 1)[0]
    row("legacy_ols_slope_ratio", "lineage / rest", "ratio", ols(L) / ols(~L), note="untransformed counts, as in the previous Fig 5")
    return s


# ----------------------------------------------------------------------------------------------------------------------
# 3. figure
# ----------------------------------------------------------------------------------------------------------------------
def draw(ax, panel, s, num):
    letter, dm, ycol, lcol, lname, rname, ylabel, hcol = panel
    L = s[lcol].values
    g, y = s.proteome_size.values.astype(float), s[ycol].values.astype(float)
    ed = edges()
    q = lambda grp, qty, an="primary": num[(num.panel == letter) & (num.analysis == an) & (num.group == grp) & (num.quantity == qty)].iloc[0]
    for m, col, mk, sz, ec, nm in ((~L, BASE_COLOR, "o", 16, "none", rname), (L, hcol, "^", 34, "black", lname)):
        pos = m & (y > 0)
        ax.scatter(g[pos], y[pos], s=sz, c=col, alpha=0.5 if mk == "o" else 0.85, linewidths=0.4 if mk == "^" else 0, edgecolors=ec, marker=mk,
                   zorder=2 if mk == "o" else 4, label=f"{nm} (n={int(m.sum())})")
        n, mg, my, se, u = bins(g[m], y[m], ed, MIN_BIN)
        ax.errorbar(mg[u], my[u], yerr=se[u], fmt="s" if mk == "o" else "D", ms=6.5, color=FIT_COLOR if mk == "o" else col, mfc="white",
                    mec=FIT_COLOR if mk == "o" else col, mew=1.4, elinewidth=1.2, capsize=0, zorder=6)
        al = q(nm, "alpha")
        lam = float(al.note.split("lambda=")[1])
        xs = np.geomspace(g[m].min(), g[m].max(), 50)
        ax.plot(xs, lam * xs ** al.estimate, "-" if mk == "^" else "--", color=col if mk == "^" else FIT_COLOR, lw=2.0 if mk == "^" else 1.6, zorder=5)
    al_l, al_r, dl, fd = q(lname, "alpha"), q(rname, "alpha"), q("lineage - rest", "delta_alpha"), q("lineage / rest", "fold_at_median_lineage_size")
    box = dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc", alpha=0.92)
    right = False                      # the fits rise to the right, so the upper left is free in every panel
    ax.text(0.965 if right else 0.035, 0.965, f"{lname}: α = {al_l.estimate:.2f} [{al_l.ci_lo:.2f}, {al_l.ci_hi:.2f}]\n{rname}: α = {al_r.estimate:.2f} [{al_r.ci_lo:.2f}, {al_r.ci_hi:.2f}]\n"
            f"Δα = {dl.estimate:+.2f} [{dl.ci_lo:+.2f}, {dl.ci_hi:+.2f}]\n×{fd.estimate:.1f} at {int(float(fd.note.split('size=')[1].split(' ')[0])):,} proteins [{fd.ci_lo:.1f}, {fd.ci_hi:.1f}]",
            transform=ax.transAxes, va="top", ha="right" if right else "left", fontsize=9.5, color=hcol if hcol != ACTINO_COLOR else "#b8621b", fontweight="bold", bbox=box)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(ylabel, fontsize=11)
    ax.grid(True, which="major", ls=":", alpha=0.35)
    ax.tick_params(labelsize=9)
    ax.set_xlim(*((1.2e3, 2.0e4) if dm == "Bacteria" else (2.0e3, 1.5e5)))   # per-domain x range: bacteria top out near 1.6e4 proteins
    ax.legend(loc="lower right", fontsize=10, frameon=True)


def main():
    global IO
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir")
    a = ap.parse_args()
    if a.dir:
        IO = Path(a.dir)
    else:
        sys.path.insert(0, str(REPO / "1_clustering"))
        from mmseqs_clusters import newest_release
        IO = Path(os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE"))) / "fig5"
    t = build_table()
    t.drop(columns=["lineage"]).to_csv(IO / "fig5_scaling_table.csv", index=False)
    rng = np.random.default_rng(SEED)
    out, binrows, subs = [], [], {}
    for p in PANELS:
        print(f"panel {p[0]}: {p[1]} {p[2]}", flush=True)
        subs[p[0]] = analyse(p, t, rng, out, binrows)
    num = pd.DataFrame(out)
    num.to_csv(IO / "fig5_scaling_numbers.csv", index=False, float_format="%.6g")
    pd.DataFrame(binrows).to_csv(IO / "fig5_scaling_bins.csv", index=False, float_format="%.6g")
    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5))
    for p, ax in zip(PANELS, axes.flat):
        draw(ax, p, subs[p[0]], num)
        ax.text(0.0, 1.02, p[0], transform=ax.transAxes, fontsize=20, fontweight="bold", ha="left", va="bottom")
    for ax in axes[1, :]:
        ax.set_xlabel("Proteome size (number of proteins)", fontsize=11)
    fig.tight_layout()
    fig.savefig(IO / "fig5_scaling_loglog.png", dpi=200, bbox_inches="tight")
    fig.savefig(IO / "fig5_scaling_loglog.pdf", bbox_inches="tight")
    print("wrote fig5_scaling_loglog.png/.pdf, fig5_scaling_numbers.csv, fig5_scaling_bins.csv, fig5_scaling_table.csv")
    pd.set_option("display.width", 220, "display.max_colwidth", 80)
    print(num[num.analysis.isin(["primary", "legacy_ols_slope_ratio"]) & num.quantity.isin(["alpha", "delta_alpha", "fold_at_median_lineage_size", "ratio"])]
          .drop(columns=["domain", "measure"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
