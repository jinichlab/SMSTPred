#!/usr/bin/env python3
"""Widen the 5 gradient rings (Protein/DNA/RNA/Metabolite/%SAMs-in-proteome) baked
into the iTOL SVG export, by rewriting their per-leaf wedge paths directly.

Geometry (verified empirically, see session notes): datasetHolder has 6
children -- 5 ring groups (1259 leaf-wedge <path> each) then 1 bar group
(1259 <g><rect> each). Every leaf's wedge in ring k sits at the SAME radius
band in mainHolder's local frame (center = local origin): ring k spans
[r_start_k, r_start_k + 25.0] with zero gap between rings (ring k+1 starts
exactly where ring k ends), confirmed by r_start/r_end having ~0 std across
all 1259 leaves. The bar dataset starts exactly at ring 4's end radius.

Each leaf-wedge path is `M P0 l (P1-P0) c c1a c1b (P2-P1) l (P3-P2) c c2a c2b
(P0-P3)`: P0/P1 at the inner/outer radius at angle a, P2/P3 at the outer/inner
radius at angle a+da (one leaf's angular slice). The l-segments are radial;
the c-segments are tangential arcs at fixed radius (c1 at r_end, c2 at
r_start). Rewriting: recompute P0..P3 from each point's ORIGINAL angle at the
NEW radius, and scale each curve's control-point offsets by (new_r/old_r) for
the radius that curve sits at -- this preserves the tangential curvature
direction while matching the new radius's larger arc.

Usage: python3 widen_rings.py [factor]   (default factor 3.0 -> 75-unit rings)
Input:  itol_tree_export_Fig3.svg, the raw iTOL export from the Zenodo record (in the data directory,
        data/ or $SMSTPRED_DATA).
Output: V4prbvTzDgPu09RO61xuyg_widened.svg next to this script (read by make_tree_figure.py).
"""
import math
import os
import re
import sys
from pathlib import Path

from lxml import etree

HERE = Path(__file__).resolve().parent
DATA_DIR = Path(os.environ.get("SMSTPRED_DATA", HERE.parents[1] / "data"))
SRC = DATA_DIR / "itol_tree_export_Fig3.svg"
N_RINGS = 5
RING_THICKNESS = 25.0          # local units, verified constant across all leaves/rings
RING0_START = 1160.596679807108  # local-frame radius where ring 0 begins (verified constant)

PATH_RE = re.compile(
    r"M(?P<p0x>-?[\d.]+),(?P<p0y>-?[\d.]+)"
    r"l(?P<l1x>-?[\d.]+),(?P<l1y>-?[\d.]+)"
    r"c(?P<c1ax>-?[\d.]+),(?P<c1ay>-?[\d.]+) (?P<c1bx>-?[\d.]+),(?P<c1by>-?[\d.]+) (?P<c1ex>-?[\d.]+),(?P<c1ey>-?[\d.]+)"
    r"l(?P<l2x>-?[\d.]+),(?P<l2y>-?[\d.]+)"
    r"c(?P<c2ax>-?[\d.]+),(?P<c2ay>-?[\d.]+) (?P<c2bx>-?[\d.]+),(?P<c2by>-?[\d.]+) (?P<c2ex>-?[\d.]+),(?P<c2ey>-?[\d.]+)$"
)


def fmt(v):
    return f"{v:.5f}"


def rewrite_path(d, old_r0, old_r1, new_r0, new_r1):
    m = PATH_RE.match(d)
    if not m:
        raise ValueError(f"unexpected path structure: {d[:120]}")
    g = {k: float(v) for k, v in m.groupdict().items()}

    p0 = (g["p0x"], g["p0y"])
    p1 = (p0[0] + g["l1x"], p0[1] + g["l1y"])
    p2 = (p1[0] + g["c1ex"], p1[1] + g["c1ey"])
    p3 = (p2[0] + g["l2x"], p2[1] + g["l2y"])

    ang_a = math.atan2(p0[1], p0[0])
    ang_b = math.atan2(p2[1], p2[0])
    ca, sa = math.cos(ang_a), math.sin(ang_a)
    cb, sb = math.cos(ang_b), math.sin(ang_b)

    np0 = (new_r0 * ca, new_r0 * sa)
    np1 = (new_r1 * ca, new_r1 * sa)
    np2 = (new_r1 * cb, new_r1 * sb)
    np3 = (new_r0 * cb, new_r0 * sb)

    k_outer = new_r1 / old_r1   # scale factor for the curve sitting at r_end
    k_inner = new_r0 / old_r0   # scale factor for the curve sitting at r_start

    new_l1 = (np1[0] - np0[0], np1[1] - np0[1])
    new_c1a = (g["c1ax"] * k_outer, g["c1ay"] * k_outer)
    new_c1b = (g["c1bx"] * k_outer, g["c1by"] * k_outer)
    new_c1e = (np2[0] - np1[0], np2[1] - np1[1])
    new_l2 = (np3[0] - np2[0], np3[1] - np2[1])
    new_c2a = (g["c2ax"] * k_inner, g["c2ay"] * k_inner)
    new_c2b = (g["c2bx"] * k_inner, g["c2by"] * k_inner)
    new_c2e = (np0[0] - np3[0], np0[1] - np3[1])

    return (
        f"M{fmt(np0[0])},{fmt(np0[1])}"
        f"l{fmt(new_l1[0])},{fmt(new_l1[1])}"
        f"c{fmt(new_c1a[0])},{fmt(new_c1a[1])} {fmt(new_c1b[0])},{fmt(new_c1b[1])} {fmt(new_c1e[0])},{fmt(new_c1e[1])}"
        f"l{fmt(new_l2[0])},{fmt(new_l2[1])}"
        f"c{fmt(new_c2a[0])},{fmt(new_c2a[1])} {fmt(new_c2b[0])},{fmt(new_c2b[1])} {fmt(new_c2e[0])},{fmt(new_c2e[1])}"
    )


def main():
    factor = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
    new_thickness = RING_THICKNESS * factor
    print(f"widening factor {factor}: ring thickness {RING_THICKNESS} -> {new_thickness} local units")

    tree = etree.parse(str(SRC))
    root = tree.getroot()
    main_el = next(e for e in root if e.get("id") == "mainHolder")
    dsh = next(e for e in main_el if e.get("id") == "datasetHolder")
    groups = list(dsh)
    assert len(groups) == N_RINGS + 1, f"expected {N_RINGS + 1} dataset groups, got {len(groups)}"

    old_bounds = [RING0_START + i * RING_THICKNESS for i in range(N_RINGS + 1)]
    new_bounds = [RING0_START + i * new_thickness for i in range(N_RINGS + 1)]
    print("old ring bounds:", [round(b, 2) for b in old_bounds])
    print("new ring bounds:", [round(b, 2) for b in new_bounds])

    for k in range(N_RINGS):
        old_r0, old_r1 = old_bounds[k], old_bounds[k + 1]
        new_r0, new_r1 = new_bounds[k], new_bounds[k + 1]
        paths = [c for c in groups[k] if etree.QName(c).localname == "path"]
        for p in paths:
            p.set("d", rewrite_path(p.get("d"), old_r0, old_r1, new_r0, new_r1))
        print(f"  ring {k}: rewrote {len(paths)} paths ({old_r0:.2f}-{old_r1:.2f} -> {new_r0:.2f}-{new_r1:.2f})")

    # Bar dataset: shift start radius out to the new ring-4 end; length (width) and
    # angular thickness (height) are data-driven and untouched.
    shift = new_bounds[N_RINGS] - old_bounds[N_RINGS]
    bar_group = groups[N_RINGS]
    n_bars = 0
    for sub_g in bar_group:
        rect = next(c for c in sub_g if etree.QName(c).localname == "rect")
        rect.set("x", fmt(float(rect.get("x")) + shift))
        n_bars += 1
    print(f"  bars: shifted {n_bars} rects outward by {shift:.2f}")

    out = HERE / "V4prbvTzDgPu09RO61xuyg_widened.svg"
    tree.write(str(out), xml_declaration=False)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
