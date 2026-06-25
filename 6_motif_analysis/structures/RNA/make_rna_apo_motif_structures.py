#!/usr/bin/env python3
"""
Apo / AlphaFold 'motif-on-the-fold' figures for RNA motifs without an RNA-bound
structure: motifs 1/4/5 on RsmF apo (3M6U), motif 8 on the RlmE AlphaFold model.
Whole protein chain A grey cartoon + motif coloured sticks (CPK heteroatoms) + label.
"""
import csv
from pathlib import Path
RNA = Path(__file__).resolve().parent
# n, source, chain, up_s, up_e, offset, hex, folder, role
MOTIFS = [
    (1, "3M6U.pdb",    "A", 172, 185, 0, "cb181d", "rna_motif1_m5Ccat", "RsmF m5C catalytic PC (apo)"),
    (4, "3M6U.pdb",    "A", 85, 101,  0, "6a51a3", "rna_motif4",        "RsmF region (apo)"),
    (5, "3M6U.pdb",    "A", 289, 297, 0, "238b45", "rna_motif5",        "RsmF region (apo)"),
    (8, "AF-Q2GDL7.pdb","A", 37, 47,  0, "e6550d", "rna_motif8_RlmE",   "RlmE basic region (AlphaFold)"),
]
def subset(src, pc):
    keep=[]
    for ln in open(RNA/src):
        rec=ln[:6]
        if rec=="CRYST1": keep.append(ln); continue
        if rec in ("ATOM  ","TER   ") and ln[21]==pc: keep.append(ln)
    keep.append("END\n"); return keep
def pml(folder,a_s,a_e,hx,role):
    return f"""# RNA {role}.  Open: pymol view.pml
load {folder}.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select motif, chain A and resi {a_s}-{a_e}
color 0x{hx}, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png {folder}.png, dpi=300
"""
def cxc(folder,a_s,a_e,hx,role,n):
    return f"""# RNA {role}.  Open: chimerax view.cxc
open {folder}.pdb
set bgColor white
hide atoms
cartoon
color /A #cccccc
color /A:{a_s}-{a_e} #{hx}
show /A:{a_s}-{a_e} atoms
style /A:{a_s}-{a_e} stick
color /A:{a_s}-{a_e} byhetero
label /A:{a_s}-{a_e} text "{{0.label_one_letter_code}}{{0.number}}" height 1.2
2dlabel text "RNA motif{n}: {role}" xpos .03 ypos .04 size 15
lighting soft
graphics silhouettes true
view
# >>> rotate, then: save {folder}.png width 2400 height 1800 supersample 3
"""
rows=[]
for n,src,pc,us,ue,off,hx,folder,role in MOTIFS:
    a_s,a_e=us+off,ue+off
    d=RNA/folder; d.mkdir(exist_ok=True)
    (d/f"{folder}.pdb").write_text("".join(subset(src,pc)))
    (d/"view.pml").write_text(pml(folder,a_s,a_e,hx,role))
    (d/"view.cxc").write_text(cxc(folder,a_s,a_e,hx,role,n))
    print(f"  motif {n} [{src} /A:{a_s}-{a_e}] -> {folder}/")
    rows.append({"motif":n,"folder":folder,"structure":src.replace('.pdb',''),
                 "type":"apo PDB" if src.startswith("3") else "AlphaFold",
                 "motif_residues":f"/A:{a_s}-{a_e}","role":role})
with open(RNA/"rna_apo_motif_structures.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("[done]")
