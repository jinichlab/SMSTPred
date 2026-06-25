#!/usr/bin/env python3
"""
Per-motif PDB + pre-styled view scripts for the RNA STREME motifs that have an
RNA-bound structure among our matches.

  2ZZM  aTrm5 (tRNA m1G37 MTase) + tRNA + SAM   -> motifs 3, 6, 7   (offset 0)
  9U48  METTL16 (U6 snRNA m6A MTase) + RNA + SAM -> motif 2          (author = UniProt-13)

Each PDB keeps: motif protein region (+/-8, author numbering) + the FULL RNA chain
+ the SAM cofactor; waters/ions/other chains stripped. view.pml / view.cxc load it
pre-styled (grey ribbon, motif coloured + sticks w/ CPK heteroatoms, RNA tan slabs,
cofactor sticks, H-bond dashes, label) -- open, rotate, save.
"""
import csv
from pathlib import Path

RNA = Path(__file__).resolve().parent
FLANK = 8

# n, pdb, prot_chain, rna_chains, cofactor, up_s, up_e, offset, hex, folder, role
MOTIFS = [
    (2, "9U48", "A", ["B"], "SAM", 148, 155, -13, "cb181d", "rna_motif2_METTL16",
     "METTL16 catalytic region (m6A snRNA)"),
    (3, "2ZZM", "A", ["B"], "SAM", 198, 212, 0, "2171b5", "rna_motif3_SAMloop",
     "Trm5 SAM-binding DxGxGxG loop"),
    (6, "2ZZM", "A", ["B"], "SAM", 45, 55, 0, "6a51a3", "rna_motif6",
     "Trm5 N-terminal region"),
    (7, "2ZZM", "A", ["B"], "SAM", 233, 240, 0, "e6550d", "rna_motif7_basic",
     "Trm5 basic RNA-binding patch (KKNIKLNK)"),
]


def subset(pdb, pc, rna_chains, cof, a_lo, a_hi):
    keep = []
    for ln in open(RNA / f"{pdb}.pdb"):
        rec = ln[:6]
        if rec == "CRYST1":
            keep.append(ln); continue
        if rec not in ("ATOM  ", "HETATM", "TER   "):
            continue
        ch = ln[21]; resname = ln[17:20].strip()
        if resname == "HOH":
            continue
        if ch == pc and rec.strip() in ("ATOM", "TER"):
            try:
                rn = int(ln[22:26])
            except ValueError:
                continue
            if a_lo <= rn <= a_hi:
                keep.append(ln)
        elif ch in rna_chains:                  # full RNA (incl. modified nts)
            keep.append(ln)
        elif rec == "HETATM" and resname == cof:  # SAM cofactor
            keep.append(ln)
    keep.append("END\n")
    return keep


def pml(folder, pc, rna_chains, a_s, a_e, hx, role):
    rnasel = "+".join(rna_chains)
    return f"""# RNA {role}.  Open: pymol view.pml
load {folder}.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
show cartoon, chain {pc}
color grey70, chain {pc}
select motif, chain {pc} and resi {a_s}-{a_e}
color 0x{hx}, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s" % (one_letter[resn], resi)
set label_size, 12
select rna, chain {rnasel}
show cartoon, rna
color wheat, rna
select cof, resn SAM
show sticks, cof
color grey60, cof
util.cnc cof
distance contacts, motif, (rna or cof), 3.6, mode=2
hide labels, contacts
color black, contacts
orient motif
deselect
# >>> rotate, then:  ray 2400,1800;  png {folder}.png, dpi=300
"""


def cxc(folder, pc, rna_chains, a_s, a_e, hx, role, n):
    rnacx = ",".join(f"/{c}" for c in rna_chains)
    return f"""# RNA {role}.  Open: chimerax view.cxc
open {folder}.pdb
set bgColor white
hide atoms
cartoon
color /{pc} #c8c8c8
color /{pc}:{a_s}-{a_e} #{hx}
show /{pc}:{a_s}-{a_e} atoms
style /{pc}:{a_s}-{a_e} stick
color /{pc}:{a_s}-{a_e} byhetero
color {rnacx} #d9b38c
nucleotides {rnacx} fill
show ligand
style ligand stick
color ligand byelement
hbonds /{pc}:{a_s}-{a_e} restrict {rnacx} | ligand reveal true color #303030 dashes 8 radius 0.05
2dlabel text "RNA motif{n}: {role}" xpos .03 ypos .04 size 15
lighting soft
graphics silhouettes true
view
# >>> rotate, then:  save {folder}.png width 2400 height 1800 supersample 3
"""


def main():
    rows = []
    for n, pdb, pc, rch, cof, us, ue, off, hx, folder, role in MOTIFS:
        a_s, a_e = us + off, ue + off
        d = RNA / folder
        d.mkdir(exist_ok=True)
        lines = subset(pdb, pc, rch, cof, a_s - FLANK, a_e + FLANK)
        (d / f"{folder}.pdb").write_text("".join(lines))
        (d / "view.pml").write_text(pml(folder, pc, rch, a_s, a_e, hx, role))
        (d / "view.cxc").write_text(cxc(folder, pc, rch, a_s, a_e, hx, role, n))
        n_at = sum(1 for x in lines if x[:4] == "ATOM")
        print(f"  motif {n} [{pdb} author {a_s}-{a_e} = UniProt {us}-{ue}] -> {folder}/  ({n_at} atoms)")
        rows.append({"motif": n, "folder": folder, "structure": pdb,
                     "protein": "aTrm5" if pdb == "2ZZM" else "METTL16",
                     "rna": "tRNA" if pdb == "2ZZM" else "snRNA",
                     "motif_residues_author": f"/{pc}:{a_s}-{a_e}",
                     "motif_residues_uniprot": f"{us}-{ue}", "role": role})
    with open(RNA / "rna_motif_structures.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("[done] RNA motif PDBs + view scripts + rna_motif_structures.csv")


if __name__ == "__main__":
    main()
