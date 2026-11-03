#!/usr/bin/env python3
"""
Motif-dark UMAP: which regions of the SAM-MTase embedding landscape are entirely
STREME-negative.

Motif status is REAL and per-protein: FIMO scanned every one of the 106,608
class members against its class's 8 STREME motifs (p<1e-4; run_fimo_all.sh ->
motif_status_all.csv). A protein is "no-motif" if it matched none of them -- no
homology propagation from the 50%-id representatives anymore.

Two panels (UMAP HDBSCAN cluster -> color):
  A  pure-no-motif UMAP clusters colored by substrate; everything else grey.
     A UMAP cluster is "pure no-motif" if every LABELED member is no-motif
     (Mixed/Unassigned ignored, lenient).
  B  ALL no-motif points colored by substrate; motif points grey
     (shows most no-motif proteins are scattered inside motif families).

Reuses the existing UMAP coordinates (no recompute).
Outputs: no_motif_umap.png + no_motif_umap_clusters.csv

  python no_motif_umap.py                         the archived figure (this folder, 5_umap/sam_mtase_umap_prostt5_coords.csv)
  python no_motif_umap.py --dir <release>/motif_analysis    a refresh release: <dir>/motif_status_all.csv (run_fimo_all.sh) on <release>/umap/sam_mtase_umap_prott5_coords.csv (5_umap/run_umap.sh),
                                                  output into <dir>
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

MA = Path(__file__).resolve().parent
COORDS = Path(__file__).resolve().parents[1] / "5_umap" / "sam_mtase_umap_prostt5_coords.csv"
MOTIF_STATUS = MA / "motif_status_all.csv"     # real per-protein FIMO status
OUT_PNG = MA / "no_motif_umap.png"

SUBSTRATE_COLORS = {"RNA": "#cb181d", "DNA": "#6a51a3", "Protein": "#238b45",
                    "Metabolite": "#e6550d", "Mixed": "#636363", "Unassigned": "#dddddd"}
GREY = "#e2e2e2"          # motif / disqualified background
SUB_ORDER = ["RNA", "DNA", "Protein", "Metabolite"]


def load_labeled():
    umap = pd.read_csv(COORDS)
    ms = pd.read_csv(MOTIF_STATUS)             # Entry, class, ..., has_motif
    status = np.where(ms["has_motif"], "motif", "no-motif")
    umap["status"] = umap["Entry"].map(dict(zip(ms["Entry"], status)))
    return umap


NEAR_DARK_FRAC = 0.90     # a UMAP cluster is "near-motif-dark" if >=90% of its
                          # labeled members are no-motif. With real per-protein FIMO
                          # NO cluster is 100% no-motif (every family has >=1 motif
                          # member) -- so we report the near-dark tail instead.


def near_motif_dark_clusters(umap, thresh=NEAR_DARK_FRAC):
    """UMAP HDBSCAN clusters whose labeled members are >= thresh no-motif."""
    keep = umap[(umap["cluster"] != -1) & umap["status"].notna()]
    g = keep.groupby("cluster")["status"]
    frac = g.apply(lambda s: (s == "no-motif").mean())
    return set(frac[frac >= thresh].index)


def scatter(ax, umap, colored_mask, letter, fg_s=3.0):
    bg = umap[~colored_mask]
    ax.scatter(bg["UMAP_1"], bg["UMAP_2"], s=1.2, c=GREY, linewidths=0, rasterized=True)
    fg = umap[colored_mask]
    cols = fg["Substrate"].map(SUBSTRATE_COLORS).fillna("#999999")
    ax.scatter(fg["UMAP_1"], fg["UMAP_2"], s=fg_s, c=cols, linewidths=0, rasterized=True)
    ax.text(0.0, 1.0, letter, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=20, fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def main():
    global COORDS, MOTIF_STATUS, OUT_PNG, MA
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=None, help="<release>/motif_analysis (default: the archived figure)")
    a = ap.parse_args()
    if a.dir:
        MA = Path(a.dir).resolve()
        COORDS = MA.parent / "umap" / "sam_mtase_umap_prott5_coords.csv"
        MOTIF_STATUS = MA / "motif_status_all.csv"
        OUT_PNG = MA / "no_motif_umap.png"
    umap = load_labeled()
    near = near_motif_dark_clusters(umap)

    in_near = umap["cluster"].isin(near)                 # panel A: near-motif-dark clusters
    is_no_motif = umap["status"] == "no-motif"           # panel B: all no-motif points

    n_near_pts = int(in_near.sum())
    n_nomotif = int(is_no_motif.sum())
    pct = int(NEAR_DARK_FRAC * 100)
    print(f"near-motif-dark clusters (>={pct}% no-motif): {len(near)} | "
          f"points in them: {n_near_pts:,} (no-motif {int((in_near & is_no_motif).sum()):,})")
    print(f"all no-motif points on UMAP: {n_nomotif:,}")

    fig, axes = plt.subplots(1, 2, figsize=(15, 7.2))
    scatter(axes[0], umap, in_near, "A", fg_s=9.0)
    scatter(axes[1], umap, is_no_motif, "B")

    handles = [mpatches.Patch(color=SUBSTRATE_COLORS[s], label=s) for s in SUB_ORDER]
    handles.append(mpatches.Patch(color=GREY, label="motif / other (grey)"))
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=10,
               frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=[0, 0.05, 1, 1.0])
    fig.savefig(OUT_PNG, dpi=200)
    print(f"[wrote] {OUT_PNG}")

    # dump the near-motif-dark cluster table for reference
    rows = []
    for c in sorted(near):
        sub = umap[umap["cluster"] == c]
        lab = sub[sub["status"].notna()]
        rows.append({"umap_cluster": c, "n_points": len(sub),
                     "n_labeled": len(lab),
                     "nomotif_frac": round((lab["status"] == "no-motif").mean(), 3),
                     "top_substrate": sub["Substrate"].mode().iat[0]})
    if rows:
        pd.DataFrame(rows).sort_values("n_points", ascending=False).to_csv(
            MA / "no_motif_umap_clusters.csv", index=False)
        print(f"[wrote] {MA/'no_motif_umap_clusters.csv'}")


if __name__ == "__main__":
    main()
