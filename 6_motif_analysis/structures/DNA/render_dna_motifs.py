#!/usr/bin/env python3
"""
Render one figure per DNA STREME motif on its DNA-bound structure (ChimeraX,
offscreen). Each figure: grey protein cartoon + tan DNA + cofactor sticks, the
motif coloured on the cartoon, and the motif residues shown as sticks with
heteroatom CPK colouring (carbon = motif colour, N = blue, O = red, S = yellow).

Output: structures/DNA/<folder>/<folder>.png  (+ the .cxc used).
"""
import subprocess
from pathlib import Path

DNA = Path(__file__).resolve().parent
CHIMERAX = "/home/nroethler/bin/chimerax"

# n, pdb, residue range, motif colour, folder, caption
MOTIFS = [
    (1, "1MHT", "/A:76-88",   "#cb181d", "motif1_PCQ",       "Motif IV catalytic PCQ (Cys81)"),
    (2, "1MHT", "/A:17-24",   "#2171b5", "motif2_SAMloop",   "Motif I SAM glycine loop"),
    (3, "1MHT", "/A:157-167", "#238b45", "motif3_TRDbasic",  "TRD DNA-contacting basic loop"),
    (4, "1MHT", "/A:274-291", "#6a51a3", "motif4_TRD",       "Target-recognition domain"),
    (5, "1MHT", "/A:300-315", "#ce1256", "motif5_TRD",       "Target-recognition domain"),
    (6, "1MHT", "/A:116-124", "#e6550d", "motif6_GluVI",     "Motif VI catalytic Glu119"),
    (7, "1G38", "/A:103-110", "#cb181d", "motif7_NPPY",      "Motif IV catalytic NPPY (Asn105)"),
    (8, "1G38", "/A:19-27",   "#2171b5", "motif8_SAMregion", "Motif I SAM region"),
]


def cxc(pdb, rng, color, png):
    drop = "delete /D,E,F\n" if pdb == "1G38" else ""   # keep one complex for TaqI
    return f"""open {DNA}/{pdb}.pdb
set bgColor white
{drop}hide atoms
cartoon
color #1 #cccccc
select nucleic
color sel #9c7a5b
~select
show ligand
style ligand stick
color ligand byelement
color {rng} {color}
show {rng} atoms
style {rng} stick
color {rng} byhetero
view {rng}
zoom 0.55
lighting soft
graphics silhouettes true
save {png} width 2000 height 1500 supersample 3
exit
"""


def main():
    for n, pdb, rng, color, folder, cap in MOTIFS:
        d = DNA / folder
        d.mkdir(exist_ok=True)
        png = d / f"{folder}.png"
        script = d / "render.cxc"
        script.write_text(cxc(pdb, rng, color, png))
        subprocess.run([CHIMERAX, "--nogui", "--offscreen", str(script)],
                       capture_output=True, timeout=300)
        print(f"  motif {n} [{pdb} {rng}] {cap}: "
              + ("OK" if png.exists() else "FAILED") + f" -> {folder}/{folder}.png")


if __name__ == "__main__":
    main()
