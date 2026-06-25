#!/usr/bin/env python3
"""
Write a pre-styled viewer script per DNA motif (PyMOL .pml and ChimeraX .cxc) that
loads <folder>.pdb and applies all the styling automatically -- motif coloured +
sticks with heteroatom CPK (N blue, O red, S yellow), grey ribbon context, DNA as
tan base-slabs, cofactor sticks, polar-contact/H-bond dashes, and a 2D label.
You just open it, rotate to taste, and save the image (the save command is left as
a comment at the bottom of each script).

Output: structures/DNA/<folder>/view.pml  and  <folder>/view.cxc
"""
from pathlib import Path

DNA = Path(__file__).resolve().parent

# n, folder, s, e, hex, cofactor, target-base, consensus, role
MOTIFS = [
    (1, "motif1_PCQ",      76, 88,  "cb181d", "SAH", "C36", "IGGPPCQGFSRAG",      "catalytic PCQ (Cys81)"),
    (2, "motif2_SAMloop",  17, 24,  "2171b5", "SAH", "C36", "PFAGSGGT",           "SAM Motif I loop"),
    (3, "motif3_TRDbasic", 157, 167,"238b45", "SAH", "C36", "YGVPQRRERVF",        "TRD DNA-contact"),
    (4, "motif4_TRD",      274, 291,"6a51a3", "SAH", "C36", "LTPREAARLQGFPDDYTF", "cofactor-proximal"),
    (5, "motif5_TRD",      300, 315,"ce1256", "SAH", "C36", "KQIGNAVPPPVARA",     "DNA-cofactor junction"),
    (6, "motif6_GluVI",    116, 124,"e6550d", "SAH", "C36", "FLLENVPGL",          "catalytic Glu119"),
    (7, "motif7_NPPY",     103, 110,"cb181d", "NEA", "6MA", "YTBPPYNT",           "catalytic NPPY (Asn105)"),
    (8, "motif8_SAMregion",19, 27,  "2171b5", "NEA", "6MA", "GZFFTPREI",          "DNA-cofactor junction"),
]


def pml(n, folder, s, e, hx, cof, tgt, cons, role):
    return f"""# {cons} - {role}.  Open:  pymol {folder}/view.pml
load {folder}.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
set valence, 1
# protein motif region as grey ribbon
show cartoon, chain A
color grey70, chain A
# the motif: colour + sticks + CPK heteroatoms
select motif, chain A and resi {s}-{e}
color 0x{hx}, motif
show sticks, motif
util.cnc motif
# DNA duplex as tan base-slabs
select dna, chain B+C
show cartoon, dna
color wheat, dna
# cofactor + target base as sticks
select cof, resn {cof}+{tgt}
show sticks, cof
color grey60, cof
util.cnc cof
# polar contacts (H-bonds) motif <-> DNA/cofactor
distance contacts, motif, (dna or cof), 3.6, mode=2
hide labels, contacts
color black, contacts
set_title m, 1, "motif{n} {cons} - {role}"
orient motif
deselect
# >>> rotate to taste, then:   ray 2400, 1800;   png motif{n}.png, dpi=300
"""


def cxc(n, folder, s, e, hx, cof, tgt, cons, role):
    rng = f"/A:{s}-{e}"
    return f"""# {cons} - {role}.  Open:  chimerax {folder}/view.cxc
open {folder}.pdb
set bgColor white
hide atoms
cartoon
color /A #b0b0b0
color {rng} #{hx}
show {rng} atoms
style {rng} stick
color {rng} byhetero
color nucleic #d9b38c
show ligand
style ligand stick
color ligand byelement
hbonds {rng} restrict nucleic | ligand reveal true color #303030 dashes 8 radius 0.05
2dlabel text "motif{n} {cons}  ({role})" xpos .03 ypos .04 size 16
lighting soft
graphics silhouettes true
view
# >>> rotate to taste, then:  save motif{n}.png width 2400 height 1800 supersample 3
"""


def main():
    for m in MOTIFS:
        n, folder = m[0], m[1]
        d = DNA / folder
        (d / "view.pml").write_text(pml(*m))
        (d / "view.cxc").write_text(cxc(*m))
        print(f"  {folder}/view.pml + view.cxc")
    print("[done] pre-styled PyMOL + ChimeraX scripts written")


if __name__ == "__main__":
    main()
