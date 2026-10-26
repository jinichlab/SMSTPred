#!/usr/bin/env python3
"""Figure 1 — study-overview pipeline (single-row flow), styled as a clean
editorial schematic: numbered stage chips + hairline rules (no saturated pills),
muted accents, thin connectors. Each of the 5 stages carries a real mini-visual.

  1 Curate     horizontal bar of the training substrate composition (imbalance)
  2 Represent  structure thumbnail -> colored embedding vector strip
  3 Train      mini cluster-disjoint MCC bar (embedding MLPs vs homology)
  4 Deploy     the tree-of-life thumbnail (Fig 3 crop)
  5 Interpret  substrate-colored UMAP + a motif-on-structure snapshot
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

BASE = Path(__file__).resolve().parent
FIG = BASE.parent
TREE = FIG / "Fig3_tree_landscape/tree_only.png"
STRUCT_REP = BASE / "biorender_components/struct_1OR8_SAM.png"   # stage 2: plain grey PRMT1 + faint SAH
STRUCT_INT = FIG / "Fig9_structures/Known/dna_motif1.png"        # stage 5: motif on structure
REPO = Path(__file__).resolve().parents[2]
COORDS = REPO / "5_umap" / "sam_mtase_umap_prostt5_coords.csv"
SUMMARY = REPO / "2_training" / "sam_classifier_comparison_output" / "sam_classifier_comparison_summary.csv"
Image.MAX_IMAGE_PIXELS = None

SUB = {"RNA": "#cb181d", "DNA": "#6a51a3", "Protein": "#238b45", "Other": "#e6550d"}
MLP_C, HOM_C = "#3b6fa0", "#c07033"       # muted embedding-MLP / homology-baseline
INK, MUTE = "#242a30", "#5c6670"
# muted, desaturated per-stage accents (chips + rules); restrained, not garish
ACCENTS = ["#3f6293", "#3d857f", "#5a8a55", "#b3833f", "#7a5c92"]
plt.rcParams.update({"font.family": "DejaVu Sans",
                     "svg.fonttype": "none"})   # keep text as editable <text> in SVG

# stage geometry (figure fraction)
XL = [0.02 + i * 0.196 for i in range(5)]
W = 0.176
XC = [x + W / 2 for x in XL]
VIZ = [[XL[i] + 0.007, 0.395, W - 0.014, 0.33] for i in range(5)]
CHIP_Y, RULE_Y, CAP_Y, ARROW_Y = 0.92, 0.865, 0.29, 0.56

TRAIN = {"DNA": 108, "RNA": 5006, "Protein": 1015, "Other": 1290}
STAGES = [
    ("1", "Curate", ["Reviewed SwissProt SAM-MTases", "4-way multi-label labels",
                     "n = 7,405 · strong class imbalance"]),
    ("2", "Represent", ["Structure-aware embeddings", "SaProt · ProstT5 · ProtT5 · ESM-IF",
                        "from AlphaFold structures · up to 1024-D"]),
    ("3", "Train + validate", ["Multi-label embedding MLP", "homology-aware cluster-disjoint CV",
                               "vs MMseqs2 / Foldseek / HMM"]),
    ("4", "Deploy", ["1,259 organisms · ~112,000 MTases", "deep ensemble → calibrated calls",
                     "+ confidence / applicability domain"]),
    ("5", "Interpret", ["Tree-of-life & taxonomic trends", "UMAP families · STREME motifs",
                        "structure mapping · discovery"]),
]


def clean_struct(path):
    """Crop the PyMOL 'License expired' watermark (red text, top-left), then trim white."""
    im = np.asarray(Image.open(path).convert("RGB"))
    corner = im[0:90, 0:520].astype(int)
    r, g, b = corner[:, :, 0], corner[:, :, 1], corner[:, :, 2]
    redrows = np.where(((r > 150) & (g < 130) & (b < 130) & (r - g > 60)).any(axis=1))[0]
    im = im[(redrows.max() + 6) if len(redrows) else 0:]
    nz = ~((im[:, :, 0] > 248) & (im[:, :, 1] > 248) & (im[:, :, 2] > 248))
    ys, xs = np.where(nz)
    return im[ys.min():ys.max(), xs.min():xs.max()]


def load_thumb(path):
    im = np.asarray(Image.open(path).convert("RGB"))
    nz = ~((im[:, :, 0] > 248) & (im[:, :, 1] > 248) & (im[:, :, 2] > 248))
    ys, xs = np.where(nz)
    return im[ys.min():ys.max(), xs.min():xs.max()]


def load_rgba(path):
    """RGBA array cropped to its non-transparent bounding box (for the plain
    grey structure render, which has a transparent background)."""
    im = np.asarray(Image.open(path).convert("RGBA"))
    ys, xs = np.where(im[:, :, 3] > 0)
    return im[ys.min():ys.max(), xs.min():xs.max()]


def _clean_axes(ax, title=None):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_edgecolor("#d5d9dd"); s.set_linewidth(0.8)
    if title:
        ax.set_title(title, fontsize=8.8, color=MUTE, pad=3)


# ---------------------------------------------------------------- viz drawers
def viz_curate(ax):
    labs = ["RNA", "Other", "Protein", "DNA"]
    vals = [TRAIN[l] for l in labs]
    y = np.arange(len(labs))[::-1]
    ax.barh(y, vals, color=[SUB[l] for l in labs], edgecolor="white", height=0.70)
    for yi, v in zip(y, vals):
        ax.text(v + 120, yi, f"{v:,}", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=9)
    ax.set_xlim(0, 6000); ax.set_xticks([])
    ax.set_title("training labels", fontsize=8.8, color=MUTE, pad=3)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_edgecolor("#c9ced3")


def viz_represent(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    try:
        st = load_rgba(STRUCT_REP)               # plain grey PRMT1 + faint SAH (RGBA)
        h, w = st.shape[:2]
        bw = 4.6; bh = bw * h / w                 # preserve aspect (landscape)
        y0 = 5.65 - bh / 2
        ax.imshow(st, extent=(0.1, 0.1 + bw, y0, y0 + bh), aspect="auto", zorder=2)
    except Exception:
        ax.add_patch(plt.Rectangle((0.4, 3.2), 3.6, 5.0, fc="#e5e5e5", ec="#bbb"))
    ax.text(2.4, 3.4, "protein structure", ha="center", va="center", fontsize=8.5, color=MUTE)
    ax.add_patch(FancyArrowPatch((4.9, 5.7), (5.5, 5.7), arrowstyle="-|>",
                                 mutation_scale=12, color=MUTE, lw=1.5))
    cmap = LinearSegmentedColormap.from_list("", ["#38548c", "#f4f4f4", "#b2382f"])
    ax.imshow(np.linspace(0, 1, 22)[None, :], extent=(5.8, 9.8, 4.8, 6.6),
              aspect="auto", cmap=cmap, zorder=2)
    ax.add_patch(plt.Rectangle((5.8, 4.8), 4.0, 1.8, fill=False, ec="#b0b4b8", lw=0.8))
    ax.text(7.8, 4.0, "embedding vector", ha="center", va="center", fontsize=8.5, color=MUTE)


def viz_train(ax):
    s = pd.read_csv(SUMMARY)
    s = s[(s.split_mode == "cluster") & (s.dataset == "full")].set_index("classifier")
    # embedding MLPs in descending MCC, then homology baselines; SaProt leads (it wins)
    picks = [("MLP_SaProt", "SaProt", MLP_C), ("MLP_ProstT5", "ProstT5", MLP_C),
             ("MLP_ProtT5", "ProtT5", MLP_C), ("MLP_ESM-IF", "ESM-IF", MLP_C),
             ("MMseqs2_1NN", "MMseqs2", HOM_C), ("Profile_HMM", "HMM", HOM_C)]
    vals = [s.loc[k, "mcc_mean"] for k, _, _ in picks]
    cols = [c for _, _, c in picks]
    x = np.arange(len(picks))
    bars = ax.bar(x, vals, color=cols, edgecolor="white", width=0.76)
    bars[1].set_edgecolor("black"); bars[1].set_linewidth(1.4)     # ProstT5 = deployed
    ax.set_ylim(0.5, 1.0); ax.set_yticks([0.5, 0.75, 1.0])
    ax.tick_params(axis="y", labelsize=7.5)
    ax.set_xticks(x); ax.set_xticklabels([n for _, n, _ in picks], rotation=32,
                                         ha="right", fontsize=6.6)
    ax.tick_params(axis="x", pad=0.5)
    ax.set_title("cluster-disjoint MCC", fontsize=8.8, color=MUTE, pad=3)
    ax.grid(axis="y", alpha=0.25)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_edgecolor("#c9ced3")


def viz_deploy(ax):
    ax.imshow(load_thumb(TREE)); ax.axis("off")
    ax.set_title("tree of life", fontsize=8.8, color=MUTE, pad=3)


def viz_interpret(ax):
    ax.axis("off")
    uax = ax.inset_axes([0.0, 0.0, 0.56, 1.0])
    sax = ax.inset_axes([0.63, 0.05, 0.37, 0.90])
    # UMAP (left)
    df = pd.read_csv(COORDS)
    if len(df) > 25000:
        df = df.sample(25000, random_state=0)
    order = ["Unassigned", "Mixed", "Other", "Protein", "DNA", "RNA"]
    palette = {**SUB, "Mixed": "#bbbbbb", "Unassigned": "#e2e2e2"}
    for sub in order:
        m = df["Substrate"] == sub
        if m.sum():
            uax.scatter(df.loc[m, "UMAP_1"], df.loc[m, "UMAP_2"],
                        s=(0.5 if sub in ("Unassigned", "Mixed") else 1.2),
                        c=palette[sub], linewidths=0, rasterized=True)
    _clean_axes(uax, "embedding UMAP")
    # motif on structure (right)
    try:
        sax.imshow(clean_struct(STRUCT_INT))
    except Exception:
        sax.add_patch(plt.Rectangle((0, 0), 1, 1, fc="#e5e5e5", ec="#bbb"))
    sax.axis("off"); sax.set_title("motif on structure", fontsize=8.8, color=MUTE, pad=3)


VIZ_FN = [viz_curate, viz_represent, viz_train, viz_deploy, viz_interpret]


def main():
    fig = plt.figure(figsize=(22, 5.0), facecolor="white")
    bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, 1); bg.set_ylim(0, 1)
    bg.axis("off"); bg.set_zorder(0)

    for i, (num, name, cap) in enumerate(STAGES):
        acc = ACCENTS[i]
        cx = XL[i] + 0.009
        # numbered chip (true circle via marker) + stage name
        bg.scatter([cx], [CHIP_Y], s=470, marker="o", color=acc,
                   edgecolors="none", zorder=4, clip_on=False)
        bg.text(cx, CHIP_Y, num, ha="center", va="center", color="white",
                fontsize=12.5, fontweight="bold", zorder=5)
        bg.text(cx + 0.015, CHIP_Y, name, ha="left", va="center", color=INK,
                fontsize=13.5, fontweight="bold")
        # hairline accent rule under the header
        bg.plot([XL[i], XL[i] + W], [RULE_Y, RULE_Y], color=acc, lw=2.0,
                solid_capstyle="butt", zorder=3)
        # viz
        VIZ_FN[i](fig.add_axes(VIZ[i]))
        # caption
        bg.text(XC[i], CAP_Y, "\n".join(cap), ha="center", va="top",
                fontsize=9.8, color="#3a4148", linespacing=1.5)
        # thin connector
        if i < 4:
            bg.add_patch(FancyArrowPatch((XL[i] + W - 0.003, ARROW_Y),
                                         (XL[i + 1] + 0.003, ARROW_Y),
                                         arrowstyle="-|>", mutation_scale=15,
                                         lw=1.5, color="#9aa3ad", zorder=4))

    out = BASE / "Fig1_flow.png"
    fig.savefig(out, dpi=200, facecolor="white")
    # vector exports: SVG keeps text/shapes editable (Illustrator/Inkscape, and
    # PowerPoint via Insert > Picture > SVG > Convert to Shape); PDF for print
    fig.savefig(BASE / "Fig1_flow.svg", facecolor="white")
    fig.savefig(BASE / "Fig1_flow.pdf", facecolor="white")
    plt.close(fig)
    print("wrote", out, "+ .svg + .pdf")


if __name__ == "__main__":
    main()
