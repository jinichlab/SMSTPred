#!/usr/bin/env python3
"""Figure S1 -- cosine-similarity basis for the 0.6 clustering threshold (Methods 2.6).
(A) All pairwise cosine similarities among the training set, before clustering.
(B) Deployment-domain max-similarity-to-training-set validation, chosen threshold marked."""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import matplotlib

BASE = os.path.dirname(os.path.abspath(__file__))
PANELS = [
    ("A", "S1A_within_training_pairwise_similarity.png"),
    ("B", "S1B_deployment_domain_validation.png"),
]

FDIR = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
F_LETTER = ImageFont.truetype(os.path.join(FDIR, "DejaVuSans-Bold.ttf"), 54)

CW, CH = 1500, 900
MARGIN, GAP = 50, 44


def trim(path):
    im = np.array(Image.open(path).convert("RGB"))
    nonwhite = ~((im[:, :, 0] > 248) & (im[:, :, 1] > 248) & (im[:, :, 2] > 248))
    ys, xs = np.where(nonwhite)
    pad = 12
    y0, y1 = max(0, ys.min() - pad), min(im.shape[0], ys.max() + pad)
    x0, x1 = max(0, xs.min() - pad), min(im.shape[1], xs.max() + pad)
    return Image.fromarray(im[y0:y1, x0:x1])


def cell(path, letter, cw, ch):
    img = trim(path)
    canvas = Image.new("RGB", (cw, ch), "white")
    s = min(cw / img.width, ch / img.height)
    nw, nh = int(img.width * s), int(img.height * s)
    img = img.resize((nw, nh), Image.LANCZOS)
    canvas.paste(img, ((cw - nw) // 2, (ch - nh) // 2))
    d = ImageDraw.Draw(canvas)
    for dx, dy in [(-2, -2), (2, -2), (-2, 2), (2, 2)]:
        d.text((10 + dx, 4 + dy), letter, font=F_LETTER, fill="white")
    d.text((10, 4), letter, font=F_LETTER, fill="black")
    return canvas


def main():
    W = 2 * CW + GAP + 2 * MARGIN
    H = CH + 2 * MARGIN
    fig = Image.new("RGB", (W, H), "white")
    for i, (letter, fname) in enumerate(PANELS):
        c = cell(os.path.join(BASE, fname), letter, CW, CH)
        x = MARGIN + i * (CW + GAP)
        y = MARGIN
        fig.paste(c, (x, y))

    out = os.path.join(BASE, "FigS1_similarity_thresholds_combined.png")
    fig.save(out, dpi=(300, 300))
    print(f"wrote {out} ({W}x{H})")


if __name__ == "__main__":
    main()
