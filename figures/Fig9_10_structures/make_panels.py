#!/usr/bin/env python3
"""Assemble each category folder's 4 structure snapshots into a labelled A/B/C/D
figure. Crops the PyMOL 'License expired' watermark (top strip) and the surrounding
whitespace, then lays the panels out 2x2 with per-panel captions."""
import argparse
import glob
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont

BASE = os.path.dirname(os.path.abspath(__file__))

# category : ordered [(filename, caption)]  (A, B, C, D)
# Snapshots use the paper palette: grey protein, class-colored motif, tan nucleic /
# wheat partner, gold small-molecule substrate, and copper (#b5651d) cofactor carbons.
FIGS = {
    "Known": [
        ("dna_motif1.png",     "M.HhaI · catalytic PC loop (DNA-1)"),
        ("rna_motif1.png",     "NSUN6 · catalytic PC loop (RNA-1)"),
        ("protein_motif5.png", "PRMT1 · THW substrate loop (Protein-5)"),
        ("metabolite_motif6.png",   "GedA O-MTase · active site (Metabolite-6)"),
    ],
    "Distal": [
        ("dna_motif4.png",     "CMT3 · junction-adjacent element (DNA-4)"),
        ("rna_motif6.png",     "NSUN6 · substrate-recognition arm (RNA-6)"),
        ("rna_motif7.png",     "aTrm5 · basic RNA-binding patch (RNA-7)"),
        ("metabolite_motif5.png",   "norcoclaurine 6-OMT · distal element (Metabolite-5)"),
    ],
}

# fonts (DejaVu ships with matplotlib)
import matplotlib
FDIR = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
F_LETTER = ImageFont.truetype(os.path.join(FDIR, "DejaVuSans-Bold.ttf"), 62)
F_CAP = ImageFont.truetype(os.path.join(FDIR, "DejaVuSans.ttf"), 30)

CW, CH = 1150, 680          # image area per cell
CAP_H = 52                  # caption bar
GAP, MARGIN = 34, 44


def clean(path):
    """crop the top watermark strip, then trim surrounding whitespace."""
    im = np.array(Image.open(path).convert("RGB"))
    # watermark = red-ish text in top-left; find its bottom row and cut below it
    corner = im[0:90, 0:520]
    r, g, b = corner[:, :, 0].astype(int), corner[:, :, 1].astype(int), corner[:, :, 2].astype(int)
    redrows = np.where(((r > 150) & (g < 130) & (b < 130) & (r - g > 60)).any(axis=1))[0]
    top = (redrows.max() + 6) if len(redrows) else 0
    im = im[top:]
    # trim whitespace: bbox of non-near-white pixels
    nonwhite = ~((im[:, :, 0] > 248) & (im[:, :, 1] > 248) & (im[:, :, 2] > 248))
    ys, xs = np.where(nonwhite)
    pad = 18
    y0, y1 = max(0, ys.min() - pad), min(im.shape[0], ys.max() + pad)
    x0, x1 = max(0, xs.min() - pad), min(im.shape[1], xs.max() + pad)
    return Image.fromarray(im[y0:y1, x0:x1])


def cell(path, letter, caption):
    img = clean(path)
    canvas = Image.new("RGB", (CW, CH + CAP_H), "white")
    # fit image into CW x CH preserving aspect
    s = min(CW / img.width, CH / img.height)
    nw, nh = int(img.width * s), int(img.height * s)
    img = img.resize((nw, nh), Image.LANCZOS)
    canvas.paste(img, ((CW - nw) // 2, (CH - nh) // 2))
    d = ImageDraw.Draw(canvas)
    # letter (top-left, with white halo for legibility over structure)
    for dx, dy in [(-2, -2), (2, -2), (-2, 2), (2, 2)]:
        d.text((14 + dx, 6 + dy), letter, font=F_LETTER, fill="white")
    d.text((14, 6), letter, font=F_LETTER, fill="black")
    # caption bar
    d.line([(0, CH), (CW, CH)], fill="#dddddd", width=1)
    font = F_CAP
    for size in range(30, 17, -1):            # shrink a long caption until it fits the cell
        font = ImageFont.truetype(os.path.join(FDIR, "DejaVuSans.ttf"), size)
        if d.textbbox((0, 0), caption, font=font)[2] <= CW - 20:
            break
    tb = d.textbbox((0, 0), caption, font=font)
    d.text(((CW - (tb[2] - tb[0])) // 2, CH + (CAP_H - (tb[3] - tb[1])) // 2 - 4),
           caption, font=font, fill="#222222")
    return canvas


def release_figs(rel):
    """The panels of a release: <release>/fig9_10/selection.csv (select_structures.py) + panels/<motif>.png (render_panels.py); returns ({category: [(png, caption)]}, folder)."""
    import sys
    sys.path.insert(0, BASE)
    from select_structures import caption
    import pandas as pd
    d = os.path.join(rel, "fig9_10")
    sel = pd.read_csv(os.path.join(d, "selection.csv"))
    figs = {cat: [(os.path.join(d, "panels", f"{r['motif']}.png"), caption(r)) for _, r in g.iterrows()] for cat, g in sel.groupby("panel_set", sort=False)}
    return figs, d


def main():
    ap = argparse.ArgumentParser(description="Fig 9 / Fig 10 layout: the archived panels (default) or a release's (--release <release>)")
    ap.add_argument("--release", default=None)
    a = ap.parse_args()
    figs, outdir = (release_figs(a.release) if a.release else (None, BASE))
    for cat, panels in (figs or FIGS).items():
        ch = CH + CAP_H
        W = 2 * CW + GAP + 2 * MARGIN
        H = 2 * ch + GAP + 2 * MARGIN
        fig = Image.new("RGB", (W, H), "white")
        for i, (fname, cap) in enumerate(panels):
            letter = chr(ord("A") + i)
            c = cell(fname if figs else os.path.join(BASE, cat, fname), letter, cap)
            col, row = i % 2, i // 2
            x = MARGIN + col * (CW + GAP)
            y = MARGIN + row * (ch + GAP)
            fig.paste(c, (x, y))
        num = 9 if cat == "Known" else 10          # Fig 9 = known sites, Fig 10 = distal sites
        out = os.path.join(outdir, f"Fig{num}_{cat}.png")
        fig.save(out, dpi=(300, 300))
        print(f"  wrote {out}  ({W}x{H})")


if __name__ == "__main__":
    main()
