#!/usr/bin/env python3
"""Fig 3 (tree landscape) drawn from the repo's own files -- no iTOL session, no SVG patching.

  python draw_tree.py [--src <folder>] [--out <folder>] [--name fig3_tree] [--gamma 0.45] [--bar-height 0.36]

--src   folder written by 4_phylo_trends/reference_tree.py: reference_tree.nwk + the itol_*.txt annotation files
        (default: <repo>/4_phylo_trends/itol_output_reference, the archived ProstT5 run; a refresh release: <release>/tree,
        written by `reference_tree.py --release`)
--out   where the figure goes (default: --src)
--gamma the radius of an inner node is (depth / max depth) ** gamma: below 1 spreads the shallow, crowded bacterial part of the
        taxonomy outwards (1 = linear, as iTOL does)

Layout (inside out): the taxonomy (Newick; unit branch lengths, so this is a TAXONOMIC tree, not a sequence phylogeny; the faint lines
join every leaf to the rim) | five rings, one leaf per angular slice: % Protein / DNA / RNA / Metabolite calls among the organism's
predicted SAM-MTases and % of the proteome that is SAM-MTase (colour scale from 0 to the largest value) | bars: number of SAM-MTases
(IPR029063) in the proteome | clade strips (Plants, Metazoa > Vertebrates, Actinomycetota; named on the figure, not in the legend). Leaves run clockwise from the top in the
order of the Newick file. Writes <name>.png (200 dpi), <name>.svg and <name>.pdf.
"""
import argparse
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

# rings inside out: (file, legend label, colour); colours come from the files themselves, these are the legend names
RINGS = [("itol_pct_protein", "% Protein-binding SAMs"), ("itol_pct_dna", "% DNA-binding SAMs"), ("itol_pct_rna", "% RNA-binding SAMs"),
         ("itol_pct_metabolite", "% Metabolite-binding SAMs"), ("itol_pct_sams", "% SAMs in proteome")]
BAR_FILE, BAR_LABEL = "itol_n_sams", "# SAM-MTases per proteome"
# clade strips: label -> (colour, track); a track is a concentric strip, Vertebrates sits outside Metazoa. Colours of the paper's legend.
CLADES = {"Plants": ("#a68d00", 0), "Metazoa": ("#0e6b74", 0), "Vertebrates": ("#7a1f4d", 1), "Actinomycetota": ("#4a2f7a", 0)}
WHITE = "#f7fbff"


# ---------------------------------------------------------------- Newick + iTOL files
class Node:
    __slots__ = ("name", "children", "parent", "depth", "angle", "radius")

    def __init__(self, parent=None):
        self.name, self.children, self.parent = "", [], parent


def parse_newick(text):
    """Nested parentheses, node labels, optional :length (lengths are ignored: every branch counts 1)."""
    s = text.strip().rstrip(";")
    root = cur = Node()
    i = 0
    label = re.compile(r"([^(),:;]*)(?::[0-9.eE+-]+)?")
    while i < len(s):
        c = s[i]
        if c == "(":
            child = Node(cur); cur.children.append(child); cur = child; i += 1
        elif c == ",":
            cur = cur.parent; child = Node(cur); cur.children.append(child); cur = child; i += 1
        elif c == ")":
            cur = cur.parent; i += 1
        else:
            m = label.match(s, i)
            cur.name = m.group(1)
            i = m.end()
    return root


def leaves_of(n, out):
    if not n.children:
        out.append(n)
    for c in n.children:
        leaves_of(c, out)
    return out


def layout(root, gamma):
    """Leaves evenly spaced clockwise from the top; an inner node sits at the middle angle of its first and last child; radius = (depth / max depth) ** gamma, leaves at 1."""
    leaves = leaves_of(root, [])
    for i, lf in enumerate(leaves):
        lf.angle = math.pi / 2 - 2 * math.pi * (i + 0.5) / len(leaves)

    def set_depth(n, d):
        n.depth = d
        for c in n.children:
            set_depth(c, d + 1)
    set_depth(root, 0)
    maxd = max(l.depth for l in leaves)

    def place(n):
        for c in n.children:
            place(c)
        if n.children:
            n.angle = (n.children[0].angle + n.children[-1].angle) / 2
    place(root)
    stack = [root]
    while stack:
        n = stack.pop()
        n.radius = 1.0 if not n.children else (n.depth / maxd) ** gamma
        stack.extend(n.children)
    return leaves


def read_dataset(path):
    """-> (header dict, {leaf: value}) of an iTOL DATASET_GRADIENT / DATASET_SIMPLEBAR file."""
    head, data, inside = {}, {}, False
    for line in open(path):
        line = line.rstrip("\n")
        if line == "DATA":
            inside = True
        elif inside and line:
            k, v = line.split("\t")[:2]
            data[k] = float(v)
        elif line and "\t" in line:
            k, v = line.split("\t", 1)
            head[k] = v
    return head, data


def read_ranges(path):
    """itol_clade_highlights.txt -> {label: (first leaf, last leaf)}"""
    out = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) == 4 and "|" in p[0]:
            out[p[3]] = tuple(p[0].split("|"))
    return out


# ---------------------------------------------------------------- drawing
def pol(r, a):
    return r * math.cos(a), r * math.sin(a)


def annulus_slice(r0, r1, a0, a1):
    return [pol(r0, a0), pol(r1, a0), pol(r1, a1), pol(r0, a1)]


def draw_tree(ax, root):
    segs, tips = [], []
    stack = [root]
    while stack:
        n = stack.pop()
        for c in n.children:
            seg = (pol(n.radius, c.angle), pol(c.radius, c.angle))
            (tips if not c.children else segs).append(seg)
            stack.append(c)
        if len(n.children) > 1:
            a0, a1 = n.children[0].angle, n.children[-1].angle
            th = np.linspace(a0, a1, max(3, int(abs(a1 - a0) / 0.02)))
            segs.append(list(zip(n.radius * np.cos(th), n.radius * np.sin(th))))
    ax.add_collection(LineCollection(tips, colors="#b0b0b0", linewidths=0.12, alpha=0.45, zorder=1))
    ax.add_collection(LineCollection(segs, colors="black", linewidths=0.5, zorder=3))


def draw_rings(ax, leaves, src, r0, width):
    """The five gradient rings; returns the outer radius and the legend data [(label, colour, vmax)]."""
    da = 2 * math.pi / len(leaves)
    r, legend = r0, []
    for stem, label in RINGS:
        head, data = read_dataset(src / f"{stem}.txt")
        vmax = max(data.values())
        cmap = LinearSegmentedColormap.from_list(stem, [WHITE, head["COLOR_MAX"]])
        polys, cols = [], []
        for lf in leaves:
            v = data.get(lf.name)
            if v is None:
                continue
            polys.append(annulus_slice(r, r + width, lf.angle - da / 2, lf.angle + da / 2))
            cols.append(cmap(v / vmax if vmax > 0 else 0))
        ax.add_collection(PolyCollection(polys, facecolors=cols, edgecolors=cols, linewidths=0.2, zorder=2))
        legend.append((label, head["COLOR_MAX"], vmax))
        r += width
    return r, legend


def draw_bars(ax, leaves, src, r0, height):
    head, data = read_dataset(src / f"{BAR_FILE}.txt")
    mx = max(data.values())
    da = 2 * math.pi / len(leaves)
    polys = [annulus_slice(r0, r0 + height * data[lf.name] / mx, lf.angle - da / 2, lf.angle + da / 2) for lf in leaves if lf.name in data]
    ax.add_collection(PolyCollection(polys, facecolors=head["COLOR"], edgecolors=head["COLOR"], linewidths=0.2, zorder=2))
    return head["COLOR"], mx


def clade_strips(ax, leaves, src, r0, strip_w):
    """One strip per clade range, a faint wedge behind the tree, and a label. Vertebrates (inside Metazoa) is on the outer track and the Metazoa label is
    put on the part of its arc that the Vertebrates do not cover, so the two names never overlap."""
    idx = {lf.name: i for i, lf in enumerate(leaves)}
    da = 2 * math.pi / len(leaves)
    ranges = {k: v for k, v in read_ranges(src / "itol_clade_highlights.txt").items() if k in CLADES}
    span = {}
    for label, (a, b) in ranges.items():
        i0, i1 = sorted((idx[a], idx[b]))
        span[label] = (i0, i1)
    for label, (i0, i1) in span.items():
        colour, track = CLADES[label]
        ang = lambda i: leaves[i].angle                                    # angles decrease with the index (clockwise)
        th = np.linspace(ang(i0) + da / 2, ang(i1) - da / 2, 80)
        ri, ro = r0 + track * (strip_w + 0.012), r0 + track * (strip_w + 0.012) + strip_w
        ax.add_patch(plt.Polygon(list(zip(ro * np.cos(th), ro * np.sin(th))) + list(zip(ri * np.cos(th[::-1]), ri * np.sin(th[::-1]))),
                                 closed=True, facecolor=colour, edgecolor="none", zorder=4))
        wedge_th = np.linspace(th[0], th[-1], 80)
        ax.add_patch(plt.Polygon([(0, 0)] + list(zip(np.cos(wedge_th), np.sin(wedge_th))), closed=True, facecolor=colour, edgecolor="none", alpha=0.10, zorder=0))
    # labels
    top = r0 + (strip_w + 0.012) + strip_w
    for label, (i0, i1) in span.items():
        colour, track = CLADES[label]
        j0, j1 = i0, i1
        if label == "Metazoa" and "Vertebrates" in span:                   # the larger part of Metazoa outside the Vertebrates
            v0, v1 = span["Vertebrates"]
            left, right = (i0, v0 - 1), (v1 + 1, i1)
            j0, j1 = max((left, right), key=lambda t: t[1] - t[0])
        mid = (leaves[j0].angle + leaves[j1].angle) / 2
        rl = top + 0.035
        x, y = pol(rl, mid)
        flip = math.cos(mid) < 0
        rot = math.degrees(mid) + (180 if flip else 0)
        ax.plot(*zip(pol(top - 0.01 if track else r0 + strip_w, mid), pol(rl - 0.005, mid)), color=colour, lw=1.0, zorder=5)
        ax.text(x, y, label, color=colour, fontsize=15, fontweight="bold", ha="right" if flip else "left", va="center", rotation=rot,
                rotation_mode="anchor", zorder=6)
    return top


def swatch(fig, x, y, w, h, colour, gradient=True):
    ax = fig.add_axes([x, y, w, h])
    if gradient:
        ax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=LinearSegmentedColormap.from_list("", [WHITE, colour]))
        ax.add_patch(Rectangle((-0.5, -0.5), 256, 1, fill=False, edgecolor="#888", lw=1.0, transform=ax.transData))
    else:
        ax.add_patch(Rectangle((0, 0), 1, 1, facecolor=colour, edgecolor="#888", lw=0.6))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def build(src, out, name, gamma, bar_h=0.36):
    root = parse_newick((src / "reference_tree.nwk").read_text())
    leaves = layout(root, gamma)
    # geometry in units of the tree radius, worked out first so that the figure height can follow the drawing
    ring_w, ring0, strip_w = 0.075, 1.03, 0.04
    r_ring = ring0 + ring_w * len(RINGS)
    top = r_ring + bar_h + 0.02 + (strip_w + 0.012) + strip_w          # outer edge of the outer clade strip
    lim = top + 0.33                                                   # half width: room for the clade names at the sides
    y_top, y_bot = lim + 0.22, -(r_ring + bar_h + 0.06)                # names stand above the 12 o'clock end; nothing hangs below the bars
    FIG_W, TREE_FRAC = 22.0, 0.62
    tree_h = FIG_W * TREE_FRAC * (y_top - y_bot) / (2 * lim)
    fig_h = tree_h / 0.96
    fig = plt.figure(figsize=(FIG_W, fig_h), facecolor="white")
    ax = fig.add_axes([0.0, 0.02, TREE_FRAC, 0.96])
    draw_tree(ax, root)
    r, legend = draw_rings(ax, leaves, src, ring0, ring_w)
    bar_colour, bar_max = draw_bars(ax, leaves, src, r, bar_h)
    clade_strips(ax, leaves, src, r + bar_h + 0.02, strip_w)
    ax.set_xlim(-lim, lim); ax.set_ylim(y_bot, y_top)
    ax.set_aspect("equal"); ax.axis("off")

    # legend, in inches from the top of the figure; the spacing is stretched towards the height of the tree (capped) and the block is centred
    LX, SW = TREE_FRAC + 0.02, 0.06
    SLBL = LX + SW + 0.014
    FS_HEAD, FS_LABEL, FS_TICK = 28, 21, 18
    SW_H, TOP0 = 0.36, 0.35
    Y = lambda inch: 1 - inch / fig_h

    def legend_pass(pad, draw, top0=TOP0):
        cur = top0

        def heading(txt):
            nonlocal cur
            if draw:
                fig.text(LX, Y(cur), txt, color="#111111", fontsize=FS_HEAD, fontweight="bold", va="top")
            cur += 0.62 + pad

        def row(txt, colour, gradient, vmax=None, note=None):
            nonlocal cur
            if draw:
                swatch(fig, LX, Y(cur + SW_H), SW, SW_H / fig_h, colour, gradient)
                fig.text(SLBL, Y(cur + SW_H / 2), txt, color="#222222", fontsize=FS_LABEL, va="center")
                if gradient:
                    fig.text(LX, Y(cur + SW_H + 0.04), "0", color="#444444", fontsize=FS_TICK, ha="left", va="top")
                    fig.text(LX + SW, Y(cur + SW_H + 0.04), f"{vmax * 100:.0f}%", color="#444444", fontsize=FS_TICK, ha="right", va="top")
                if note:
                    fig.text(SLBL, Y(cur + SW_H + 0.02), note, color="#555555", fontsize=FS_TICK, va="top")
            cur += SW_H + (0.46 if gradient else 0.16 + (0.30 if note else 0)) + pad

        heading("Rings (inner to outer)")
        for label, colour, vmax in legend:
            row(label, colour, True, vmax)
        cur += 0.10
        heading("Outer bars")
        row(BAR_LABEL, bar_colour, False, note=f"longest bar = {bar_max:.0f}")
        return cur

    n_slots = 2 + len(legend) + 1
    target = 0.955 * fig_h - 0.1
    pad = min(0.18, max(0.0, (target - legend_pass(0.0, False)) / n_slots))      # spread out to the tree's height, but never gappy
    height = legend_pass(pad, False) - TOP0
    legend_pass(pad, True, top0=max(TOP0, (fig_h - height) / 2))                  # and centred beside the tree
    n = len(leaves)

    out.mkdir(parents=True, exist_ok=True)
    for ext, kw in (("png", {"dpi": 200}), ("svg", {}), ("pdf", {})):
        fig.savefig(out / f"{name}.{ext}", facecolor="white", **kw)
    plt.close(fig)
    print(f"{n} leaves; ring maxima " + ", ".join(f"{l.split()[0]}{l.split()[1] if l.startswith('% S') else ''} {v:.3f}" for l, _, v in legend) + f"; bar max {bar_max:.0f}")
    print(f"wrote {out}/{name}.png, .svg, .pdf")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=str(REPO / "4_phylo_trends" / "itol_output_reference"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--name", default="fig3_tree")
    ap.add_argument("--gamma", type=float, default=0.45)
    ap.add_argument("--bar-height", type=float, default=0.36, help="height of the longest SAM-count bar, in units of the tree radius")
    a = ap.parse_args()
    build(Path(a.src), Path(a.out or a.src), a.name, a.gamma, a.bar_height)


if __name__ == "__main__":
    main()
