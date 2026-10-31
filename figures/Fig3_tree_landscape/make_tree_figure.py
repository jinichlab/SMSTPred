#!/usr/bin/env python3
"""Build the Fig 3 tree panel from the iTOL SVG export (white background).

The SVG is the good source: tree-topology branches are stroke #000000 (they were
invisible on the black PNG export -- only the grey leaf connectors showed), so on
a white background the real phylogeny appears in black. We (1) strip the iTOL
legend + scale groups, (2) rasterize on white with cairosvg, (3) crop to content,
(4) compose with a compact custom legend keyed to the paper palette, and
(5) optionally bracket taxonomic clades (see HIGHLIGHTS).

CLADE HIGHLIGHTS
----------------
Each leaf label in the SVG carries its angular position as transform="rotate(a)"
and radius as x; the tree center + orientation come from the mainHolder transform
(translate . rotate . scale). We push each leaf anchor through the full affine
chain into cropped-raster pixel coordinates, match leaf label -> canon(Organism)
-> taxonomic lineage, group the matching leaves into contiguous angular arcs, and
draw an iTOL-style clade strip just outside the bar band with a label + a faint
accent wedge. Add a lineage by name in HIGHLIGHTS; a polyphyletic pick simply
yields multiple strips."""
from pathlib import Path
import copy
import io
import os
import re
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle, Wedge
from PIL import Image
from lxml import etree
import cairosvg

Image.MAX_IMAGE_PIXELS = None
BASE = Path(__file__).resolve().parent
# Gradient rings widened 3x from the raw iTOL export via widen_rings.py (the
# native export's rings are only 25 local units thick -- hard to tell apart).
SVG = BASE / "V4prbvTzDgPu09RO61xuyg_widened.svg"
SVG_RAW = Path(os.environ.get("SMSTPRED_DATA", BASE.parents[1] / "data")) / "itol_tree_export_Fig3.svg"  # raw iTOL export (Zenodo)
RENDER_W = 7000                                  # rasterize width (px) of the original 1512-unit canvas
PX_PER_UNIT = RENDER_W / 1512.0                  # raster scale (px per SVG unit) kept for the enlarged canvas
PROBE_MARGIN = 1500                              # SVG units added on every side to find the tree's true extent
DROP_GROUPS = {"legendHolder", "scaleHolder", "timeScaleHolder", "wrapperHolder"}

WHITE = "#f7fbff"          # iTOL gradient COLOR_MIN (low end)
# Ordered inner -> outer to match the concentric ring stacking on the tree
# (radii measured from this SVG export: Protein 129, DNA 132, RNA 135, Metabolite 137,
# proteome 140). "% SAMs in proteome" is the outermost ring, then the outer bars.
RINGS = [("% Protein-binding SAMs", "#238b45"),
         ("% DNA-binding SAMs", "#6a51a3"),
         ("% RNA-binding SAMs", "#cb181d"),
         ("% Metabolite (small-molecule) SAMs", "#e6550d"),
         ("% SAMs in proteome", "#2171b5")]
BAR_COLOR = "#2b8cbe"

# Clade highlights: (lineage substring, display label, color). Substring is matched
# against the organism's full "Taxonomic lineage"; monophyletic picks draw one strip,
# polyphyletic picks draw several. Add entries freely, e.g.
#   ("Vertebrata", "Vertebrates", "#1f78b4"),
#   ("Actinomycetota", "Actinomycetes", "#e6550d"),
ORG_CSV = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[2] / "data")) / "selected_1259_organisms.csv"
# The clade ranges are now baked into the iTOL SVG export (native TREE_COLORS
# ranges from reference_tree.py's itol_clade_highlights.txt), so the matplotlib
# overlay is off to avoid drawing them twice. Re-populate to add non-taxonomic
# highlights the SVG doesn't already carry.
HIGHLIGHTS = []

# Legend key for the clade ranges baked into the SVG. (taxon in "Taxonomic
# lineage", label, hex color); colors MUST match the ranges in the export
# (itol_clade_highlights.txt). Entries are re-ordered for the legend by each
# clade's angular position in the tree, so the key reads in tree order.
CLADE_KEY = [
    ("Viridiplantae",  "Plants",        "#a68d00"),
    ("Metazoa",        "Metazoa",        "#0e6b74"),
    ("Vertebrata",     "Vertebrates",   "#7a1f4d"),
    ("Actinomycetota", "Actinomycetota", "#4a2f7a"),
]


# ------------------------------------------------------------------ name matching
def canon(name):
    name = re.sub(r"\s+", "_", str(name).strip())
    name = re.sub(r"[^A-Za-z0-9_\-]+", "", name)
    return re.sub(r"_+", "_", name).strip("_") or "unnamed"


# ---- SVG affine helpers (3x3, row-vector convention p @ M) ----
def _T(tx, ty): return np.array([[1, 0, 0], [0, 1, 0], [tx, ty, 1]], float)
def _S(s):      return np.array([[s, 0, 0], [0, s, 0], [0, 0, 1]], float)
def _R(deg):
    a = math.radians(deg); c, s = math.cos(a), math.sin(a)
    return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]], float)


def _parse_transform(s):
    M = np.eye(3)
    for fn, args in re.findall(r"(translate|rotate|scale)\(([^)]*)\)", s or ""):
        v = [float(x) for x in re.split(r"[,\s]+", args.strip()) if x]
        if fn == "translate":
            M = _T(v[0], v[1] if len(v) > 1 else 0) @ M
        elif fn == "scale":
            M = _S(v[0]) @ M
        elif fn == "rotate":
            M = _R(v[0]) @ M
    return M


def contiguous_arcs(angles, gap_deg=2.5):
    """Group sorted angles (deg, 0..360) into arcs whose neighbors are < gap apart,
    merging across the 0/360 wrap. Returns list of (theta1, theta2) in degrees."""
    if not angles:
        return []
    a = sorted(angles)
    arcs, start, prev = [], a[0], a[0]
    for ang in a[1:]:
        if ang - prev > gap_deg:
            arcs.append((start, prev)); start = ang
        prev = ang
    arcs.append((start, prev))
    if len(arcs) > 1 and (arcs[0][0] + 360 - arcs[-1][1]) < gap_deg:
        arcs[0] = (arcs[-1][0], arcs[0][1] + 360); arcs.pop()
    return arcs


def _raster(root, vbox, width):
    """Rasterize the (already stripped) SVG on white over the user-unit window vbox."""
    r = copy.deepcopy(root)
    r.set("viewBox", " ".join(f"{v:.4f}" for v in vbox))
    r.set("width", f"{vbox[2]:.4f}"); r.set("height", f"{vbox[3]:.4f}")
    png = cairosvg.svg2png(bytestring=etree.tostring(r), output_width=int(width),
                           background_color="white")
    return Image.open(io.BytesIO(png)).convert("RGB")


def render_svg():
    """Strip legend/scale groups, rasterize on white, crop to content, and recover
    the tree geometry (center, leaf tip radius, per-leaf angle, bar-band radius) in
    cropped-raster pixel coordinates for clade highlighting.

    The iTOL canvas (1512 x 764 units) is smaller than the tree: the longest outer
    bars run past its edges and used to be clipped flat. We first probe a much
    larger canvas to find the tree's true extent, then rasterize exactly that region
    at the original scale."""
    tree = etree.parse(str(SVG)); root = tree.getroot()
    vb = [float(x) for x in re.split(r"[,\s]+", root.get("viewBox"))]
    main = next((_parse_transform(e.get("transform")) for e in root.iter()
                 if e.get("id") == "mainHolder"), np.eye(3))
    # leaf labels: <text> with a rotate transform + numeric x, inside mainHolder
    leaves = []
    for e in root.iter():
        if etree.QName(e).localname != "text":
            continue
        tr = e.get("transform")
        if not tr or "rotate" not in tr or e.get("x") is None or not (e.text or "").strip():
            continue
        p = np.array([float(e.get("x")), float(e.get("y") or 0), 1.0]) @ (_parse_transform(tr) @ main)
        leaves.append((e.text.strip(), p[0], p[1]))
    center_vb = np.array([0, 0, 1.0]) @ main

    for el in [el for el in root.iter() if el.get("id") in DROP_GROUPS]:
        el.getparent().remove(el)
    print(f"  dropped groups: {DROP_GROUPS & {el.get('id') for el in tree.iter()} or 'all'}")
    probe_box = [vb[0] - PROBE_MARGIN, vb[1] - PROBE_MARGIN,
                 vb[2] + 2 * PROBE_MARGIN, vb[3] + 2 * PROBE_MARGIN]
    pw = 1800
    pa = np.asarray(_raster(root, probe_box, pw))
    pys, pxs = np.where(~((pa[:, :, 0] > 250) & (pa[:, :, 1] > 250) & (pa[:, :, 2] > 250)))
    ps = pw / probe_box[2]; pad_u = 12.0
    vb_old = vb
    vb = [probe_box[0] + pxs.min() / ps - pad_u, probe_box[1] + pys.min() / ps - pad_u,
          (pxs.max() - pxs.min()) / ps + 2 * pad_u, (pys.max() - pys.min()) / ps + 2 * pad_u]
    print(f"  canvas {vb_old[2]:.0f}x{vb_old[3]:.0f} units -> tree extent {vb[2]:.0f}x{vb[3]:.0f} units "
          f"(x {vb[0]:.0f}..{vb[0] + vb[2]:.0f}, y {vb[1]:.0f}..{vb[1] + vb[3]:.0f})")
    im = _raster(root, vb, round(vb[2] * PX_PER_UNIT))
    scale = im.size[0] / vb[2]
    a = np.asarray(im)
    content = ~((a[:, :, 0] > 244) & (a[:, :, 1] > 244) & (a[:, :, 2] > 244))
    ys, xs = np.where(content); pad = 24
    box = (max(0, xs.min() - pad), max(0, ys.min() - pad),
           min(a.shape[1], xs.max() + pad), min(a.shape[0], ys.max() + pad))
    crop = im.crop(box)
    crop.save(BASE / "tree_only.png")
    print(f"  rasterized {im.size} -> tree crop {crop.size}")

    def to_px(vx, vy):
        return (vx - vb[0]) * scale - box[0], (vy - vb[1]) * scale - box[1]
    cx, cy = to_px(*center_vb[:2])
    leaf_ang = []
    leaf_r = None
    for nm, vx, vy in leaves:
        px, py = to_px(vx, vy)
        leaf_ang.append((nm, math.degrees(math.atan2(py - cy, px - cx)) % 360))
        leaf_r = math.hypot(px - cx, py - cy)          # all tips share one radius
    # bar-band radius: outside ~80% of per-angle bar tips (rare long bars cross it)
    ca = np.asarray(crop)
    cont = ~((ca[:, :, 0] > 244) & (ca[:, :, 1] > 244) & (ca[:, :, 2] > 244))
    yy, xx = np.where(cont)
    ang = np.degrees(np.arctan2(yy - cy, xx - cx)) % 360
    rad = np.hypot(xx - cx, yy - cy)
    bins = np.clip(np.floor(ang).astype(int), 0, 359)
    per_bin_max = np.array([rad[bins == b].max() if (bins == b).any() else np.nan
                            for b in range(360)])
    band_r = float(np.nanpercentile(per_bin_max, 80))
    geom = dict(center=(cx, cy), leaf_r=leaf_r, band_r=band_r, leaf_ang=leaf_ang)
    return crop, geom


def build_highlights(geom):
    """Resolve each HIGHLIGHTS entry to (label, color, arcs) using the lineage CSV."""
    if not HIGHLIGHTS:
        return []
    df = pd.read_csv(ORG_CSV)
    out = []
    for taxon, label, color in HIGHLIGHTS:
        ids = {canon(o) for o, lin in zip(df["Organism"], df["Taxonomic lineage"])
               if isinstance(lin, str) and taxon in lin}
        angs = [ang for nm, ang in geom["leaf_ang"] if canon(nm) in ids]
        arcs = contiguous_arcs(angs)
        print(f"  highlight {taxon}: {len(angs)} leaves -> {len(arcs)} arc(s)")
        out.append((label, color, arcs))
    return out


def clade_order(geom):
    """Return [(label, color, n_leaves), ...] for CLADE_KEY, ordered by each clade's
    mean angular position in the tree (ascending pixel-frame angle = the clades read
    left-to-right across the top), so the legend key matches the tree."""
    if not CLADE_KEY:
        return []
    df = pd.read_csv(ORG_CSV)
    items = []
    for taxon, label, color in CLADE_KEY:
        ids = {canon(o) for o, lin in zip(df["Organism"], df["Taxonomic lineage"])
               if isinstance(lin, str) and taxon in lin}
        angs = [math.radians(a) for nm, a in geom["leaf_ang"] if canon(nm) in ids]
        if not angs:
            print(f"  clade key: '{taxon}' has no leaves in tree - skipping")
            continue
        mean = math.degrees(math.atan2(sum(math.sin(a) for a in angs) / len(angs),
                                       sum(math.cos(a) for a in angs) / len(angs))) % 360
        items.append((mean, label, color, len(angs)))
    items.sort(key=lambda t: t[0])
    return [(label, color, n) for _, label, color, n in items]


def draw_highlights(ax, geom, highlights):
    """Clade strip just outside the bar band + faint accent wedge + label."""
    cx, cy = geom["center"]; leaf_r = geom["leaf_r"]
    strip_in = geom["band_r"] + leaf_r * 0.03
    strip_w = leaf_r * 0.05
    strip_out = strip_in + strip_w
    for label, color, arcs in highlights:
        for t1, t2 in arcs:
            ax.add_patch(Wedge((cx, cy), leaf_r * 1.01, t1, t2, width=leaf_r * 0.46,
                               facecolor=color, alpha=0.10, edgecolor="none", zorder=4))
            ax.add_patch(Wedge((cx, cy), strip_out, t1, t2, width=strip_w,
                               facecolor=color, edgecolor="none", zorder=6))
            mid = math.radians((t1 + t2) / 2.0)
            dx, dy = math.cos(mid), math.sin(mid)
            lead = leaf_r * 0.14
            lx0, ly0 = cx + strip_out * dx, cy + strip_out * dy
            lx1, ly1 = cx + (strip_out + lead) * dx, cy + (strip_out + lead) * dy
            ax.plot([lx0, lx1], [ly0, ly1], color=color, lw=1.6, zorder=6, clip_on=False)
            ha = "right" if dx < 0 else "left"
            ax.text(lx1 + (-1 if dx < 0 else 1) * leaf_r * 0.02, ly1, label,
                    color=color, fontsize=13.5, fontweight="bold", ha=ha, va="center",
                    linespacing=1.15, zorder=7, clip_on=False)


def swatch(fig, x, y, w, h, color, gradient=True):
    ax = fig.add_axes([x, y, w, h])
    if gradient:
        cmap = LinearSegmentedColormap.from_list("", [WHITE, color])
        ax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=cmap)
        ax.add_patch(Rectangle((-0.5, -0.5), 256, 1, fill=False, edgecolor="#888", lw=1.0,
                               transform=ax.transData))
    else:
        ax.add_patch(Rectangle((0, 0), 1, 1, facecolor=color, edgecolor="#888", lw=0.6))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


FIG_W = 22.0          # inches; the manuscript scales this to page width (~0.3x), so legend text is sized
TREE_FRAC = 0.60      # share of the figure width given to the tree
FS_HEAD, FS_LABEL, FS_TICK = 36, 31, 23   # legend fonts (pt at 22 in wide; ~10.6 / 9.1 / 6.8 pt at 6.5 in)


def compose(tree, geom):
    W, H = tree.size
    fig_h = FIG_W * TREE_FRAC * H / W / 0.96
    fig = plt.figure(figsize=(FIG_W, fig_h), facecolor="white")
    axt = fig.add_axes([0.003, 0.02, TREE_FRAC, 0.96])
    axt.imshow(tree); axt.axis("off")
    draw_highlights(axt, geom, build_highlights(geom))

    # No title / organism-info block here (Figure 3 caption in the manuscript text already
    # covers both). Legend layout is in inches from the top so the large text never collides.
    LX, SW = TREE_FRAC + 0.015, 0.068            # legend left edge, swatch width (figure fraction)
    SLBL = LX + SW + 0.014                       # label x
    SWH = 0.58                                   # swatch height (in)
    clades = clade_order(geom)
    labels = [n for n, _ in RINGS] + ["# SAM-MTases per proteome"] + [lb for lb, _, _ in clades]

    # largest label font (<= FS_LABEL) at which the longest label still fits inside the figure
    renderer = fig.canvas.get_renderer()
    avail_pt = (0.985 - SLBL) * FIG_W * 72
    fs_label = FS_LABEL
    while fs_label > 16:
        t = fig.text(0, 0, max(labels, key=len), fontsize=fs_label)
        w_pt = t.get_window_extent(renderer).width * 72 / fig.dpi
        t.remove()
        if w_pt <= avail_pt:
            break
        fs_label -= 1
    print(f"  legend fonts: heading {FS_HEAD}, label {fs_label}, tick {FS_TICK} pt; figure {FIG_W:.1f} x {fig_h:.1f} in")

    Y = lambda inch: 1 - inch / fig_h           # inches from top -> figure fraction
    cur = 0.35

    def heading(txt):
        nonlocal cur
        fig.text(LX, Y(cur), txt, color="#111111", fontsize=FS_HEAD, fontweight="bold", va="top")
        cur += 0.80

    def row(name, color, gradient):
        nonlocal cur
        swatch(fig, LX, Y(cur + SWH), SW, SWH / fig_h, color, gradient=gradient)
        fig.text(SLBL, Y(cur + SWH / 2), name, color="#222222", fontsize=fs_label, va="center")
        if gradient:
            fig.text(LX, Y(cur + SWH + 0.04), "0", color="#444444", fontsize=FS_TICK, ha="left", va="top")
            fig.text(LX + SW, Y(cur + SWH + 0.04), "max", color="#444444", fontsize=FS_TICK, ha="right", va="top")
            cur += SWH + 0.46 + 0.22
        else:
            cur += SWH + 0.32

    heading("Gradient rings")
    for name, color in RINGS:
        row(name, color, True)
    cur += 0.15
    heading("Outer bars")
    row("# SAM-MTases per proteome", BAR_COLOR, False)
    cur += 0.20
    if clades:
        heading("Highlighted clades")
        for label, color, _n in clades:
            row(label, color, False)
    print(f"  legend height used {cur:.1f} of {fig_h:.1f} in")

    out = BASE / "fig3_tree_landscape.png"
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    crop, geom = render_svg()
    compose(crop, geom)
