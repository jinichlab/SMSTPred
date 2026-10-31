#!/usr/bin/env python3
"""Figure 7 — enzyme-family clusters underlying the tree-of-life trends.
Row 1 (A/B) = characterized families recovered de novo and anchored to
reviewed SwissProt members (PRMT/CARM; DNMT3A/3B). Row 2 (C/D) = uncharacterized
'dark-matter' families with no reviewed members (plant and actinomycete
small-molecule enzymes); the four outlined clusters in panel C (c110, c551, c667,
c714) are the ones examined in Figure 8 (CLEAN / Horizyn reaction profiles)."""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import matplotlib

BASE = os.path.dirname(os.path.abspath(__file__))
PANELS = [
    ("A", "A_metazoa_protein_clusters.png"),
    ("B", "B_vertebrate_dna_clusters.png"),
    ("C", "C_plant_metabolite_dark_clusters.png"),
    ("D", "D_actino_metabolite_dark_clusters.png"),
]

FDIR = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
F_LETTER = ImageFont.truetype(os.path.join(FDIR, "DejaVuSans-Bold.ttf"), 72)

CW, CH = 1740, 1180        # A-D cell content box
MARGIN, GAP, TITLE_H = 50, 44, 20


def trim(path):
    im = np.array(Image.open(path).convert("RGB"))
    nonwhite = ~((im[:, :, 0] > 248) & (im[:, :, 1] > 248) & (im[:, :, 2] > 248))
    ys, xs = np.where(nonwhite)
    pad = 12
    y0, y1 = max(0, ys.min() - pad), min(im.shape[0], ys.max() + pad)
    x0, x1 = max(0, xs.min() - pad), min(im.shape[1], xs.max() + pad)
    return Image.fromarray(im[y0:y1, x0:x1])


LETTER_BAND = 90  # reserved strip at the top of every cell, so the big panel
                   # letter never lands on top of the panel's own content


def cell(path, letter, cw, ch):
    img = trim(path)
    content_h = ch - LETTER_BAND
    canvas = Image.new("RGB", (cw, ch), "white")
    s = min(cw / img.width, content_h / img.height)
    nw, nh = int(img.width * s), int(img.height * s)
    img = img.resize((nw, nh), Image.LANCZOS)
    canvas.paste(img, ((cw - nw) // 2, LETTER_BAND + (content_h - nh) // 2))
    d = ImageDraw.Draw(canvas)
    for dx, dy in [(-2, -2), (2, -2), (-2, 2), (2, 2)]:
        d.text((10 + dx, 4 + dy), letter, font=F_LETTER, fill="white")
    d.text((10, 4), letter, font=F_LETTER, fill="black")
    return canvas


def main():
    W = 2 * CW + GAP + 2 * MARGIN
    H = TITLE_H + 2 * CH + GAP + 2 * MARGIN
    fig = Image.new("RGB", (W, H), "white")
    for i, (letter, fname) in enumerate(PANELS):
        c = cell(os.path.join(BASE, fname), letter, CW, CH)
        col, row = i % 2, i // 2
        x = MARGIN + col * (CW + GAP)
        y = MARGIN + TITLE_H + row * (CH + GAP)
        fig.paste(c, (x, y))

    out = os.path.join(BASE, "Fig7_clusters_combined.png")
    fig.save(out, dpi=(300, 300))
    print(f"wrote {out} ({W}x{H})")


if __name__ == "__main__":
    main()
