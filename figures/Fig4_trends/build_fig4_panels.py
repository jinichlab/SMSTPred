#!/usr/bin/env python3
"""Figure 4 source panels (A-D) -- substrate-load trends across the tree of life,
rebuilt with NO embedded titles (axis labels + in-plot stat annotations only;
the manuscript caption carries the description). Replaces the four PNGs that
`make_fig4_trends.py` composites into Fig4_trends_combined.png.

Panel A: Protein-substrate load, Metazoa vs unicellular vs non-metazoan
         multicellular eukaryotes (reproduces repo/4_phylo_trends/protein_cellularity_violin.py).
Panel B: DNA- and Protein-substrate load, Vertebrata vs invertebrate Metazoa.
Panel C: Other-substrate load, Plants (Viridiplantae) vs non-plant eukaryotes.
Panel D: Other-substrate load, Actinomycetota vs other Bacteria.

Source: Stage-3 ProstT5 deployment predictions (repo/3_predictions/predictions_output/
*_embeddings_predictions.csv, multi-label call_* columns) joined to
data/selected_1259_organisms.csv
(Taxonomic lineage, is_metazoa, is_vertebrate).

Outputs (this dir): A_protein_3group_violin.png, B_vertebrate_invertebrate_violins.png,
C_plant_other_violin.png, D_actino_other_violin.png
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

REPO = Path(__file__).resolve().parents[2]
PRED_DIR = REPO / "3_predictions" / "predictions_output"
META_CSV = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[2] / "data")) / "selected_1259_organisms.csv"
PRED_SUFFIX = "_embeddings_predictions.csv"
OUT_DIR = Path(__file__).resolve().parent
MIN_N = 5
YEASTS = ("Saccharomycotina", "Schizosaccharomyces", "Taphrinomycotina")


def organism_to_stem(name):
    return str(name).replace(" ", "_").replace("/", "_")


def cellularity(lineage):
    if not isinstance(lineage, str) or "Eukaryota" not in lineage:
        return None
    if "Metazoa" in lineage:
        return "Metazoa"
    if "Embryophyta" in lineage:
        return "Non-metazoan multicellular"
    if "Phaeophyceae" in lineage or "Rhodophyta" in lineage:
        return "Non-metazoan multicellular"
    if "Oomyc" in lineage or "Peronospor" in lineage:
        return "Non-metazoan multicellular"
    if "Fungi" in lineage:
        if any(y in lineage for y in YEASTS):
            return "Unicellular eukaryote"
        return "Non-metazoan multicellular"
    return "Unicellular eukaryote"


def plant_group(lineage):
    if not isinstance(lineage, str) or "Eukaryota" not in lineage:
        return None
    return "Plants" if "Viridiplantae" in lineage else "Non-plant eukaryotes"


def bacteria_group(lineage):
    if not isinstance(lineage, str) or "Bacteria" not in lineage:
        return None
    return "Actinomycetota" if "Actinomycetota" in lineage else "Other bacteria"


def stars(p):
    return ("****" if p < 1e-4 else "***" if p < 1e-3 else
            "**" if p < 0.01 else "*" if p < 0.05 else "ns")


def build_org_table(meta):
    stem_to_lineage = {organism_to_stem(o): lin for o, lin in
                        zip(meta["Organism"], meta["Taxonomic lineage"])}
    stem_to_org = {organism_to_stem(o): o for o in meta["Organism"]}
    stem_to_vert = {organism_to_stem(o): v for o, v in zip(meta["Organism"], meta["is_vertebrate"])}
    stem_to_metazoa = {organism_to_stem(o): v for o, v in zip(meta["Organism"], meta["is_metazoa"])}

    rows = []
    for f in sorted(PRED_DIR.glob(f"*{PRED_SUFFIX}")):
        try:
            d = pd.read_csv(f, low_memory=False)
        except Exception:
            continue
        if d.empty:
            continue
        stem = f.name[:-len(PRED_SUFFIX)]
        n_total = len(d)
        row = {"Organism": stem_to_org.get(stem, stem), "lineage": stem_to_lineage.get(stem),
               "is_vertebrate": stem_to_vert.get(stem, False), "is_metazoa": stem_to_metazoa.get(stem, False),
               "n_total": n_total}
        for cls in ["Protein", "DNA", "Other"]:
            col = f"call_{cls}"
            n = int((pd.to_numeric(d[col], errors="coerce") == 1).sum()) if col in d.columns else 0
            row[cls] = n
            row[f"pct_{cls}"] = n / n_total * 100 if n_total else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def violin(ax, data_by_group, groups, palette, ylab):
    arrays = [data_by_group[g] for g in groups]
    parts = ax.violinplot(arrays, positions=range(len(groups)), showmedians=True, showextrema=False, widths=0.85)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(palette[groups[i]]); body.set_edgecolor("black"); body.set_alpha(0.75)
    parts["cmedians"].set_color("black")
    rng = np.random.default_rng(0)
    for i in range(len(groups)):
        x = i + rng.uniform(-0.10, 0.10, size=len(arrays[i]))
        ax.scatter(x, arrays[i], s=8, color="black", alpha=0.25, linewidths=0, zorder=3)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g.replace(" ", "\n") for g in groups], fontsize=24)
    ax.set_ylabel(ylab, fontsize=24); ax.grid(axis="y", linestyle=":", alpha=0.3)
    ax.tick_params(axis="both", labelsize=21)
    vals = np.concatenate(arrays) if any(len(a) for a in arrays) else np.array([0.0])
    ymax = max(float(vals.max()) * 1.15, 5.0)
    ax.set_ylim(-ymax * 0.10, ymax * 1.15)
    for i, g in enumerate(groups):
        ax.text(i, -0.22, f"n={len(data_by_group[g])}", ha="center", va="top", fontsize=19,
                color="#555", transform=ax.get_xaxis_transform())
    return ymax


def pairwise_brackets(ax, data_by_group, groups, ymax, kruskal_box=True):
    pairs = []
    for (i, gi), (j, gj) in combinations(enumerate(groups), 2):
        a, b = data_by_group[gi], data_by_group[gj]
        if len(a) >= 3 and len(b) >= 3:
            U, p = mannwhitneyu(a, b, alternative="two-sided")
            pairs.append((i, j, p))
    y0, dy = ymax * 1.08, ymax * 0.16
    for k, (i, j, p) in enumerate(pairs):
        y = y0 + k * dy
        ax.plot([i, i, j, j], [y, y + ymax * 0.018, y + ymax * 0.018, y], color="black", lw=1.6)
        ax.text((i + j) / 2, y + ymax * 0.026, stars(p), ha="center", va="bottom",
                fontsize=26, fontweight="bold")
    top = y0 + max(len(pairs), 1) * dy
    if len(groups) > 2 and kruskal_box:
        H, p_kw = kruskal(*[data_by_group[g] for g in groups])
        kruskal_y = y0 + len(pairs) * dy + ymax * 0.05
        ax.text((len(groups) - 1) / 2, kruskal_y, f"Kruskal H={H:.1f}, p={p_kw:.1e}",
                ha="center", va="bottom", fontsize=19, color="#444",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#ccc", alpha=0.9))
        top = kruskal_y + ymax * 0.16
    ax.set_ylim(top=top)
    return pairs


def panel_A(org):
    groups = ["Metazoa", "Unicellular eukaryote", "Non-metazoan multicellular"]
    palette = {"Metazoa": "#d95f02", "Unicellular eukaryote": "#7570b3", "Non-metazoan multicellular": "#1b9e77"}
    org = org.copy()
    org["grp"] = org["lineage"].apply(cellularity)
    sub = org.dropna(subset=["grp"])
    sub = sub[sub["n_total"] >= MIN_N]
    print("Panel A counts:\n", sub["grp"].value_counts().to_string())

    fig, axes = plt.subplots(1, 2, figsize=(16, 8.5), constrained_layout=True)
    data_count = {g: sub.loc[sub["grp"] == g, "Protein"].to_numpy() for g in groups}
    yc = violin(axes[0], data_count, groups, palette, "Protein SAM-MTases per organism")
    pairwise_brackets(axes[0], data_count, groups, yc)
    data_pct = {g: sub.loc[sub["grp"] == g, "pct_Protein"].to_numpy() for g in groups}
    yp = violin(axes[1], data_pct, groups, palette, "% of organism's SAM-MTases")
    pairwise_brackets(axes[1], data_pct, groups, yp)
    fig.savefig(OUT_DIR / "A_protein_3group_violin.png", dpi=200)
    plt.close(fig)
    print("[wrote] A_protein_3group_violin.png")


def panel_B(org):
    groups = ["Vertebrata", "Invertebrate\nMetazoa"]
    org = org.copy()
    met = org[org["is_metazoa"] == True].copy()
    met["grp"] = met["is_vertebrate"].map({True: "Vertebrata", False: "Invertebrate\nMetazoa"})
    met = met[met["n_total"] >= MIN_N]
    print("Panel B counts:\n", met["grp"].value_counts().to_string())
    palette = {"Vertebrata": "#b2604a", "Invertebrate\nMetazoa": "#e8b98a"}

    fig, axes = plt.subplots(2, 2, figsize=(16, 15.5), constrained_layout=True)
    for row, cls in enumerate(["DNA", "Protein"]):
        for col, metric in enumerate(["count", "pct"]):
            ax = axes[row, col]
            key = cls if metric == "count" else f"pct_{cls}"
            data = {g: met.loc[met["grp"] == g, key].to_numpy() for g in groups}
            ylab = f"{cls} SAM-MTases" if metric == "count" else "% of organism's SAM-MTases"
            ymax = violin(ax, data, groups, palette, ylab)
            _, p = mannwhitneyu(data[groups[0]], data[groups[1]])
            ax.plot([0, 0, 1, 1], [ymax * 1.07, ymax * 1.09, ymax * 1.09, ymax * 1.07], color="black", lw=1.6)
            ax.text(0.5, ymax * 1.12, f"{stars(p)}  (p={p:.1e})", ha="center", va="bottom", fontsize=26, fontweight="bold")
            ax.set_ylim(top=ymax * 1.34)
    fig.savefig(OUT_DIR / "B_vertebrate_invertebrate_violins.png", dpi=200)
    plt.close(fig)
    print("[wrote] B_vertebrate_invertebrate_violins.png")


def panel_C(org):
    groups = ["Plants", "Non-plant\neukaryotes"]
    palette = {"Plants": "#2ca25f", "Non-plant\neukaryotes": "#9e9ac8"}
    org = org.copy()
    org["grp"] = org["lineage"].apply(plant_group).replace({"Non-plant eukaryotes": "Non-plant\neukaryotes"})
    sub = org.dropna(subset=["grp"])
    sub = sub[sub["n_total"] >= MIN_N]
    print("Panel C counts:\n", sub["grp"].value_counts().to_string())

    fig, axes = plt.subplots(1, 2, figsize=(16, 13), constrained_layout=True)
    for ax, key, ylab in [(axes[0], "Other", "Other-substrate SAM-MTases per organism"),
                           (axes[1], "pct_Other", "% of organism's SAM-MTases that are 'Other'")]:
        data = {g: sub.loc[sub["grp"] == g, key].to_numpy() for g in groups}
        ymax = violin(ax, data, groups, palette, ylab)
        _, p = mannwhitneyu(data[groups[0]], data[groups[1]])
        ax.plot([0, 0, 1, 1], [ymax * 1.07, ymax * 1.09, ymax * 1.09, ymax * 1.07], color="black", lw=1.6)
        ax.text(0.5, ymax * 1.12, f"{stars(p)}  (p={p:.1e})", ha="center", va="bottom", fontsize=26, fontweight="bold")
        ax.set_ylim(top=ymax * 1.34)
    fig.savefig(OUT_DIR / "C_plant_other_violin.png", dpi=200)
    plt.close(fig)
    print("[wrote] C_plant_other_violin.png")


def panel_D(org):
    groups = ["Actinomycetota", "Other\nbacteria"]
    palette = {"Actinomycetota": "#d95f02", "Other\nbacteria": "#bdbdbd"}
    org = org.copy()
    org["grp"] = org["lineage"].apply(bacteria_group).replace({"Other bacteria": "Other\nbacteria"})
    sub = org.dropna(subset=["grp"])
    sub = sub[sub["n_total"] >= MIN_N]
    print("Panel D counts:\n", sub["grp"].value_counts().to_string())

    fig, axes = plt.subplots(1, 2, figsize=(16, 13), constrained_layout=True)
    for ax, key, ylab in [(axes[0], "Other", "Other-substrate SAM-MTases per organism"),
                           (axes[1], "pct_Other", "% of organism's SAM-MTases that are 'Other'")]:
        data = {g: sub.loc[sub["grp"] == g, key].to_numpy() for g in groups}
        ymax = violin(ax, data, groups, palette, ylab)
        _, p = mannwhitneyu(data[groups[0]], data[groups[1]])
        ax.plot([0, 0, 1, 1], [ymax * 1.07, ymax * 1.09, ymax * 1.09, ymax * 1.07], color="black", lw=1.6)
        ax.text(0.5, ymax * 1.12, f"{stars(p)}  (p={p:.1e})", ha="center", va="bottom", fontsize=26, fontweight="bold")
        ax.set_ylim(top=ymax * 1.34)
    fig.savefig(OUT_DIR / "D_actino_other_violin.png", dpi=200)
    plt.close(fig)
    print("[wrote] D_actino_other_violin.png")


def main():
    meta = pd.read_csv(META_CSV)
    org = build_org_table(meta)
    panel_A(org)
    panel_B(org)
    panel_C(org)
    panel_D(org)


if __name__ == "__main__":
    main()
