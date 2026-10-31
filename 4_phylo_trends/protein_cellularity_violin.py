#!/usr/bin/env python3
"""
Per-organism Protein-substrate SAM-MTase load across eukaryotic cellularity
groups, on the ProstT5 deployment predictions.

Unit = one eukaryotic organism (>= MIN_N SAM-MTases). Each organism is assigned
a cellularity group from its UniProt lineage (multicellularity is polyphyletic,
so the call is rule-based, not a single rank):

  - Metazoa                      multicellular animals
  - Unicellular eukaryote        protists, algae, unicellular yeasts
  - Non-metazoan multicellular   land plants, filamentous fungi, brown/red
                                 algae, oomycetes (plant-dominated)

Substrate calls come from the Stage-3 ProstT5 deployment predictions
(../3_predictions/predictions_output/{stem}_embeddings_predictions.csv, the
multi-label call_Protein column). Two panels: Protein-MTase count and % per
organism, with Kruskal-Wallis + pairwise Mann-Whitney significance.

Headline: Metazoa carry significantly more Protein-substrate MTases (count and
%) than either unicellular or non-metazoan multicellular eukaryotes.

Outputs (this directory): protein_3group_violin.png, protein_cellularity_stats.csv
Env: smstpred (`make env`).
"""
import os
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

PRED_DIR = Path(__file__).resolve().parent.parent / "3_predictions" / "predictions_output"
META_CSV = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[1] / "data")) / "selected_1259_organisms.csv"
PRED_SUFFIX = "_embeddings_predictions.csv"
OUT_DIR = Path(__file__).resolve().parent
OUT_PNG = OUT_DIR / "protein_3group_violin.png"
OUT_STATS = OUT_DIR / "protein_cellularity_stats.csv"

SUBSTRATE = "Protein"
MIN_N = 5
GROUPS = ["Metazoa", "Unicellular eukaryote", "Non-metazoan multicellular"]
YEASTS = ("Saccharomycotina", "Schizosaccharomyces", "Taphrinomycotina")
PALETTE = {
    "Metazoa": "#ed8c44",
    "Unicellular eukaryote": "#4578aa",
    "Non-metazoan multicellular": "#618d43",
}


def organism_to_stem(name):
    return str(name).replace(" ", "_").replace("/", "_")


def cellularity(lineage):
    """Rule-based uni/multicellular call for eukaryotes; None for non-eukaryotes."""
    if not isinstance(lineage, str) or "Eukaryota" not in lineage:
        return None
    if "Metazoa" in lineage:
        return "Metazoa"
    if "Embryophyta" in lineage:                               # land plants
        return "Non-metazoan multicellular"
    if "Phaeophyceae" in lineage or "Rhodophyta" in lineage:   # brown / red algae
        return "Non-metazoan multicellular"
    if "Oomyc" in lineage or "Peronospor" in lineage:          # water molds
        return "Non-metazoan multicellular"
    if "Fungi" in lineage:
        if any(y in lineage for y in YEASTS):                  # unicellular yeasts
            return "Unicellular eukaryote"
        return "Non-metazoan multicellular"                    # filamentous fungi
    return "Unicellular eukaryote"                             # remaining microeuks


def stars(p):
    return ("****" if p < 1e-4 else "***" if p < 1e-3 else
            "**" if p < 0.01 else "*" if p < 0.05 else "ns")


def build_org_table():
    meta = pd.read_csv(META_CSV)
    stem_to_lineage = {organism_to_stem(o): lin for o, lin in
                       zip(meta["Organism"], meta["Taxonomic lineage"])}
    stem_to_org = {organism_to_stem(o): o for o in meta["Organism"]}

    rows = []
    for f in sorted(PRED_DIR.glob(f"*{PRED_SUFFIX}")):
        try:
            d = pd.read_csv(f, low_memory=False)
        except Exception as e:
            print(f"  skip {f.name}: {e}", flush=True)
            continue
        if d.empty or "call_Protein" not in d.columns:
            continue
        stem = f.name[:-len(PRED_SUFFIX)]
        n_prot = int((pd.to_numeric(d["call_Protein"], errors="coerce") == 1).sum())
        rows.append({"Organism": stem_to_org.get(stem, stem),
                     "lineage": stem_to_lineage.get(stem),
                     "n_total": len(d), SUBSTRATE: n_prot})

    org = pd.DataFrame(rows)
    org["coarse3"] = org["lineage"].apply(cellularity)
    org[f"pct_{SUBSTRATE}"] = org[SUBSTRATE] / org["n_total"] * 100
    org = org[org["n_total"] >= MIN_N].dropna(subset=["coarse3"]).copy()
    print(f"eukaryotic organisms (>= {MIN_N} SAM-MTases): {len(org):,}", flush=True)
    print(org["coarse3"].value_counts().to_string(), flush=True)
    return org


def auto_ymax(data_by_group, min_value=5.0, pad=1.15):
    vals = np.concatenate([v for v in data_by_group.values() if len(v)])
    return max(float(np.max(vals)) * pad, min_value) if len(vals) else min_value


def violin(ax, data_by_group, groups, title, ylab, ymax):
    arrays = [data_by_group[g] for g in groups]
    parts = ax.violinplot(arrays, positions=range(len(groups)),
                          showmedians=True, showextrema=False, widths=0.85)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(PALETTE[groups[i]]); body.set_edgecolor("black"); body.set_alpha(0.7)
    parts["cmedians"].set_color("black")
    rng = np.random.default_rng(0)
    for i in range(len(groups)):
        x = i + rng.uniform(-0.10, 0.10, size=len(arrays[i]))
        ax.scatter(x, arrays[i], s=6, color="black", alpha=0.25, linewidths=0, zorder=3)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g.replace(" ", "\n") for g in groups], fontsize=9)
    ax.set_ylim(-ymax * 0.04, ymax * 1.20)
    ax.set_ylabel(ylab); ax.set_title(title); ax.grid(axis="y", linestyle=":", alpha=0.3)
    for i, g in enumerate(groups):
        ax.text(i, -0.10, f"n={len(data_by_group[g])}", ha="center", va="top",
                fontsize=7, color="#555", transform=ax.get_xaxis_transform())


def brackets_and_stats(ax, data_by_group, groups, ymax, metric, rows):
    H, p_kw = kruskal(*[data_by_group[g] for g in groups])
    ax.text(0.02, 0.97, f"Kruskal H={H:.1f}\np={p_kw:.1e}", transform=ax.transAxes,
            va="top", ha="left", fontsize=7, color="#555",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
    pairs = []
    for (i, gi), (j, gj) in combinations(enumerate(groups), 2):
        a, b = data_by_group[gi], data_by_group[gj]
        U, p = (mannwhitneyu(a, b, alternative="two-sided") if len(a) >= 3 and len(b) >= 3
                else (np.nan, np.nan))
        pairs.append((i, j, p))
        rows.append({"metric": metric, "group_a": gi, "group_b": gj,
                     "median_a": float(np.median(a)) if len(a) else np.nan,
                     "median_b": float(np.median(b)) if len(b) else np.nan,
                     "n_a": len(a), "n_b": len(b),
                     "U": float(U) if not np.isnan(U) else np.nan,
                     "p_mw": float(p) if not np.isnan(p) else np.nan,
                     "kruskal_H": float(H), "kruskal_p": float(p_kw)})
    y0, dy = ymax * 1.03, ymax * 0.06
    for k, (i, j, p) in enumerate([pp for pp in pairs if not np.isnan(pp[2])]):
        y = y0 + k * dy
        ax.plot([i, i, j, j], [y, y + ymax * 0.01, y + ymax * 0.01, y], color="black", lw=1.0)
        ax.text((i + j) / 2, y + ymax * 0.015, stars(p), ha="center", va="bottom", fontsize=9)


def main():
    org = build_org_table()
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(max(8, 2.7 * len(GROUPS)), 4.8))

    data_count = {g: org.loc[org["coarse3"] == g, SUBSTRATE].to_numpy() for g in GROUPS}
    yc = auto_ymax(data_count, 5.0)
    violin(axes[0], data_count, GROUPS, f"{SUBSTRATE}  --  count per organism",
           f"{SUBSTRATE} SAM-MTases per organism", yc)
    brackets_and_stats(axes[0], data_count, GROUPS, yc, "count", rows)

    data_pct = {g: org.loc[org["coarse3"] == g, f"pct_{SUBSTRATE}"].to_numpy() for g in GROUPS}
    yp = auto_ymax(data_pct, 10.0)
    violin(axes[1], data_pct, GROUPS, f"{SUBSTRATE}  --  % per organism",
           "% of organism's SAM-MTases", yp)
    brackets_and_stats(axes[1], data_pct, GROUPS, yp, "pct", rows)

    fig.suptitle("Protein-substrate SAM-MTases per organism (ProstT5):\n"
                 "Metazoa vs unicellular vs non-metazoan multicellular eukaryotes",
                 y=1.00, fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(OUT_STATS, index=False)
    print(f"wrote {OUT_PNG.name} and {OUT_STATS.name}", flush=True)


if __name__ == "__main__":
    main()
