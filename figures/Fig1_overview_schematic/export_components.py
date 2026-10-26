#!/usr/bin/env python3
"""Export each Figure-1 sub-panel as a standalone, transparent-background PNG for
BioRender. You upload these individual assets to BioRender and rebuild the frame
(numbered chips, accent rules, connector arrows, stage headers, captions) natively
with BioRender shapes/text -- BioRender cannot import the composite as editable
elements, so the sub-panels are provided atomically instead.

Reuses the drawers/constants in build_fig1_flow.py so the pieces match the figure.
Data-driven panels (bar charts, UMAP) are rendered onto a transparent canvas;
raster thumbnails (structures, tree) are keyed white->alpha so they composite
cleanly over any BioRender background.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
import build_fig1_flow as F     # reuse viz_curate/viz_train, palette, paths

OUT = F.BASE / "biorender_components"
OUT.mkdir(exist_ok=True)
DPI = 300


def save_panel(fn, name, w, h):
    """Render a drawer(ax) onto a tight, transparent canvas."""
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([0.16, 0.16, 0.80, 0.74])
    fn(ax)
    fig.savefig(OUT / name, dpi=DPI, transparent=True, bbox_inches="tight",
                pad_inches=0.04)
    plt.close(fig)
    print("  ", name)


def white_to_alpha(rgb, thr=244):
    rgb = np.asarray(rgb).astype(np.uint8)
    a = np.where((rgb[:, :, 0] > thr) & (rgb[:, :, 1] > thr) & (rgb[:, :, 2] > thr),
                 0, 255).astype(np.uint8)
    return np.dstack([rgb, a])


def save_image(arr, name):
    Image.fromarray(white_to_alpha(arr)).save(OUT / name)
    print("  ", name)


def save_embedding_strip():
    fig = plt.figure(figsize=(3.2, 1.0)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    cmap = LinearSegmentedColormap.from_list("", ["#38548c", "#f4f4f4", "#b2382f"])
    ax.imshow(np.linspace(0, 1, 22)[None, :], aspect="auto", cmap=cmap)
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, fill=False,
                               ec="#b0b4b8", lw=1.2))
    fig.savefig(OUT / "02b_embedding_strip.png", dpi=DPI, transparent=True,
                bbox_inches="tight", pad_inches=0.02)
    plt.close(fig); print("   02b_embedding_strip.png")


def save_umap():
    df = pd.read_csv(F.COORDS)
    if len(df) > 40000:
        df = df.sample(40000, random_state=0)
    order = ["Unassigned", "Mixed", "Other", "Protein", "DNA", "RNA"]
    palette = {**F.SUB, "Mixed": "#bbbbbb", "Unassigned": "#e2e2e2"}
    fig = plt.figure(figsize=(3.2, 3.2)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    for sub in order:
        m = df["Substrate"] == sub
        if m.sum():
            ax.scatter(df.loc[m, "UMAP_1"], df.loc[m, "UMAP_2"],
                       s=(1.2 if sub in ("Unassigned", "Mixed") else 3.0),
                       c=palette[sub], linewidths=0, rasterized=True)
    fig.savefig(OUT / "05a_umap.png", dpi=DPI, transparent=True,
                bbox_inches="tight", pad_inches=0.02)
    plt.close(fig); print("   05a_umap.png")


def main():
    print(f"[export] -> {OUT}")
    # data-driven vector-ish panels (keep their text/labels)
    save_panel(F.viz_curate, "01_curate_bars.png", 4.2, 2.6)
    save_panel(F.viz_train, "03_mcc_bars.png", 4.2, 3.0)
    # atomic raster assets (white keyed to transparent)
    save_embedding_strip()
    save_umap()
    try:
        save_image(F.load_thumb(F.TREE), "04_tree.png")
    except Exception as e:
        print("   [skip] tree:", e)
    for path, name in [(F.STRUCT_REP, "02a_structure.png"),
                       (F.STRUCT_INT, "05b_motif_structure.png")]:
        try:
            save_image(F.clean_struct(path), name)
        except Exception as e:
            print(f"   [skip] {name}:", e)
    print("[done] upload these to BioRender; rebuild chips/arrows/headers natively.")


if __name__ == "__main__":
    main()
