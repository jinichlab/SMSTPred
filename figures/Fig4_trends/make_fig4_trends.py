#!/usr/bin/env python3
"""Figure 4 — substrate-load trends across the tree of life (four violin panels,
2x2, labelled A/B/C/D). Row 1 = animal expansions (protein/DNA); row 2 = plant
and actinomycete small-molecule ('Metabolite') expansions. These are the statistical
trends that the UMAP (Figure 6) and family clusters (Figure 7) then dissect."""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import matplotlib

BASE = os.path.dirname(os.path.abspath(__file__))
PANELS = [
    ("A", "A_protein_3group_violin.png"),
    ("B", "B_vertebrate_invertebrate_violins.png"),
    ("C", "C_plant_metabolite_violin.png"),
    ("D", "D_actino_metabolite_violin.png"),
]
FDIR = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
F_LETTER = ImageFont.truetype(os.path.join(FDIR, "DejaVuSans-Bold.ttf"), 72)

CW, CH = 1740, 1180        # cell content box
MARGIN, GAP, TITLE_H = 50, 44, 20


def trim(path):
    im = np.array(Image.open(path).convert("RGB"))
    nonwhite = ~((im[:, :, 0] > 248) & (im[:, :, 1] > 248) & (im[:, :, 2] > 248))
    ys, xs = np.where(nonwhite)
    pad = 12
    y0, y1 = max(0, ys.min() - pad), min(im.shape[0], ys.max() + pad)
    x0, x1 = max(0, xs.min() - pad), min(im.shape[1], xs.max() + pad)
    return Image.fromarray(im[y0:y1, x0:x1])


def cell(path, letter):
    img = trim(path)
    canvas = Image.new("RGB", (CW, CH), "white")
    s = min(CW / img.width, CH / img.height)
    nw, nh = int(img.width * s), int(img.height * s)
    img = img.resize((nw, nh), Image.LANCZOS)
    canvas.paste(img, ((CW - nw) // 2, (CH - nh) // 2))
    d = ImageDraw.Draw(canvas)
    for dx, dy in [(-2, -2), (2, -2), (-2, 2), (2, 2)]:
        d.text((10 + dx, 4 + dy), letter, font=F_LETTER, fill="white")
    d.text((10, 4), letter, font=F_LETTER, fill="black")
    return canvas


def main():
    W = 2 * CW + GAP + 2 * MARGIN
    H = TITLE_H + 2 * CH + GAP + 2 * MARGIN
    fig = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(fig)
    for i, (letter, fname) in enumerate(PANELS):
        c = cell(os.path.join(BASE, fname), letter)
        col, row = i % 2, i // 2
        x = MARGIN + col * (CW + GAP)
        y = MARGIN + TITLE_H + row * (CH + GAP)
        fig.paste(c, (x, y))
    out = os.path.join(BASE, "Fig4_trends_combined.png")
    fig.save(out, dpi=(300, 300))
    print(f"wrote {out} ({W}x{H})")


if __name__ == "__main__":
    main()
