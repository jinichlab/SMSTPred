#!/usr/bin/env python3
"""
Interaction-focused figure per DNA STREME motif.

  - Motif region shown as a thin ribbon (+/-8 flanking residues, grey; motif span
    coloured) so the secondary-structure context is legible.
  - Motif residues' atoms as sticks, heteroatoms CPK (N blue, O red, S yellow).
  - Binding partner only:
      DNA-binding / TRD motifs -> contacting DNA (slabs for context within 7 A;
        the directly-contacted nucleotides, <3.8 A, shown as base STICKS so base
        interactions are visible).
      SAM-binding motifs       -> cofactor (SAH / AdoMet analog), sticks.
      catalytic motifs         -> both.
  - H-bonds motif<->partner drawn as dashes.
Per-motif `orient` lets us nudge the camera so the interaction faces the viewer.

Output: structures/DNA/<folder>/<folder>_interaction.png
"""
import subprocess
from pathlib import Path

DNA = Path(__file__).resolve().parent
CHIMERAX = "/home/nroethler/bin/chimerax"
FLANK = 8

# n, pdb, chain, start, end, colour, folder, show_dna, show_lig, orient
MOTIFS = [
    (1, "1MHT", "A", 76, 88,  "#cb181d", "motif1_PCQ",       True,  True,  ""),
    (2, "1MHT", "A", 17, 24,  "#2171b5", "motif2_SAMloop",   False, True,  ""),
    (3, "1MHT", "A", 157, 167,"#238b45", "motif3_TRDbasic",  True,  False, ""),
    (4, "1MHT", "A", 274, 291,"#6a51a3", "motif4_cofactor",  False, True,  ""),
    (5, "1MHT", "A", 300, 315,"#ce1256", "motif5_junction",  True,  True,  ""),
    (6, "1MHT", "A", 116, 124,"#e6550d", "motif6_GluVI",     True,  False, ""),
    (7, "1G38", "A", 103, 110,"#cb181d", "motif7_NPPY",      True,  True,  ""),
    (8, "1G38", "A", 19, 27,  "#2171b5", "motif8_junction",  True,  True,  ""),
]


def cxc(pdb, chain, s, e, color, png, show_dna, show_lig, orient):
    rng = f"/{chain}:{s}-{e}"
    ribbon = f"/{chain}:{max(1, s-FLANK)}-{e+FLANK}"
    near = f"(({rng}) :< 7) & nucleic"
    contact = f"(({rng}) :< 3.8) & nucleic"
    drop = ["delete /D,E,F"] if pdb == "1G38" else []
    out = [f"open {DNA}/{pdb}.pdb", "set bgColor white", *drop,
           "hide atoms", "hide cartoon",
           # motif region ribbon for context
           f"show {ribbon} cartoon", f"color {ribbon} #c8c8c8", f"color {rng} {color}",
           # motif atoms as sticks, CPK heteroatoms
           f"show {rng} atoms", f"style {rng} stick", f"color {rng} byhetero"]
    partners = []
    if show_dna:
        out += [f"show {near} cartoon", f"color {near} #9c7a5b",        # base slabs (context)
                f"show {contact} atoms", f"style {contact} stick",       # contacted base -> sticks
                f"nucleotides {contact} atoms",
                f"color {contact} #9c7a5b", f"color {contact} byhetero"]
        partners.append(near)
    if show_lig:
        out += ["show ligand atoms", "style ligand stick", "color ligand byelement"]
        partners.append("ligand")
    partner = " | ".join(f"({p})" for p in partners)
    out += [f"hbonds {rng} restrict {partner} reveal true color #303030 radius 0.05 dashes 8",
            f"view {rng} | {partner}", "zoom 0.8"]
    if orient:
        out.append(orient)
    out += ["lighting soft", "graphics silhouettes true",
            f"save {png} width 2000 height 1500 supersample 3", "exit"]
    return "\n".join(out)


def main():
    for n, pdb, ch, s, e, color, folder, sd, sl, orient in MOTIFS:
        d = DNA / folder
        d.mkdir(exist_ok=True)
        png = d / f"{folder}_interaction.png"
        script = d / "render_interaction.cxc"
        script.write_text(cxc(pdb, ch, s, e, color, png, sd, sl, orient))
        subprocess.run([CHIMERAX, "--nogui", "--offscreen", str(script)],
                       capture_output=True, timeout=300)
        print(f"  motif {n} [{folder}]: " + ("OK" if png.exists() else "FAILED"))


if __name__ == "__main__":
    main()
