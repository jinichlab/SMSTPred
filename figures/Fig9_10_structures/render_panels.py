#!/usr/bin/env python3
"""Render the structure snapshots of Figs 9 and 10 headlessly with ChimeraX.

  python render_panels.py --release <release> [--chimerax ~/bin/chimerax]

Reads <release>/fig9_10/selection.csv (select_structures.py) and, for each selected motif, the view.cxc that motif_structure_pick.py wrote next to its subset PDB
(<release>/motif_analysis/structure_search/structures/<Class>/<folder>/). The styling is that script's (grey enzyme, class-coloured motif, tan nucleic acid / wheat peptide,
copper SAM carbons, hydrogen bonds of the motif); here the in-picture caption is dropped (the panel caption is drawn by make_panels.py) and the camera is set on the motif and the ligands
instead of on the whole enzyme. Writes <release>/fig9_10/panels/<motif>.png (2000 x 1500).
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


def render_script(view_cxc, png, chain, rng, size=(2000, 1500)):
    keep = []
    for ln in Path(view_cxc).read_text().splitlines():
        if not ln.strip() or ln.startswith("#") or ln.startswith("2dlabel"):
            continue
        keep.append(ln.strip())
    keep = [("view /%s:%s | ligand pad 0.6" % (chain, rng)) if ln == "view" else ln for ln in keep]
    keep += [f"save {png} width {size[0]} height {size[1]} supersample 3", "exit"]
    return "\n".join(keep) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", required=True)
    ap.add_argument("--chimerax", default=os.environ.get("CHIMERAX", str(Path.home() / "bin" / "chimerax")))
    a = ap.parse_args()
    rel = Path(a.release)
    sel = pd.read_csv(rel / "fig9_10" / "selection.csv")
    base = rel / "motif_analysis" / "structure_search" / "structures"
    pdir = rel / "fig9_10" / "panels"
    pdir.mkdir(exist_ok=True)
    for _, r in sel.iterrows():
        folder = base / r["folder"]
        png = pdir / f"{r['motif']}.png"
        script = pdir / f"_{r['motif']}.cxc"
        script.write_text(f"cd {folder}\n" + render_script(folder / "view.cxc", png, r["chain"], r["range"]))
        p = subprocess.run([a.chimerax, "--offscreen", "--nogui", "--exit", "--script", str(script)], capture_output=True, text=True, timeout=600)
        if not png.exists():
            sys.exit(f"ChimeraX made no image for {r['motif']}:\n{p.stdout[-800:]}\n{p.stderr[-800:]}")
        script.unlink()
        print(f"  {r['motif']:13} -> {png}", flush=True)


if __name__ == "__main__":
    main()
