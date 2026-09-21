#!/usr/bin/env python3
"""IG of Q9V2G1 (no-motif RNA-MTase) mapped onto its tRNA-bound structure 5WT1.
IG written to B-factor (chain A); keep protein A + tRNA (C,F) + SAH. view scripts:
grey protein, top-10 IG residues as red sticks + labels, RNA tan slabs, SAH sticks."""
from pathlib import Path
import pandas as pd
H=Path(__file__).resolve().parent
ig=pd.read_csv(H.parent.parent/"attribution"/"ig_q9v2g1"/"Q9V2G1__RNA_residue_attr.csv")
igmap=dict(zip(ig.position.astype(int),ig.ig.astype(float)))
top=ig.sort_values("ig",ascending=False).head(10)
topres=list(top.position.astype(int)); topcsv=",".join(map(str,topres))
KEEP={"A","C","F"}; DROPHET={"HOH","GOL"}
out=[]
for ln in open(H/"5WT1.pdb"):
    rec=ln[:6]
    if rec=="CRYST1": out.append(ln); continue
    if rec in ("ATOM  ","HETATM","TER   "):
        ch=ln[21]; rn3=ln[17:20].strip()
        if ch not in KEEP: continue
        if rec=="HETATM" and rn3 in DROPHET: continue
        if rec=="ATOM  " and ch=="A":
            try: r=int(ln[22:26])
            except: out.append(ln); continue
            b=igmap.get(r,0.0)*100
            ln=ln[:60]+f"{b:6.2f}"+ln[66:]
        elif rec in ("ATOM  ","HETATM"):
            ln=ln[:60]+f"{0.0:6.2f}"+ln[66:]
        out.append(ln)
out.append("END\n")
(H/"Q9V2G1_IG.pdb").write_text("".join(out))
(H/"view.cxc").write_text(f"""# Q9V2G1 (no-motif RNA-MTase) IG on tRNA-bound 5WT1.  Open: chimerax view.cxc
open Q9V2G1_IG.pdb
set bgColor white
hide atoms
cartoon
color /A #cccccc
color /A:{topcsv} #cb181d
show /A:{topcsv} atoms
style /A:{topcsv} stick
color /A:{topcsv} byhetero
label /A:{topcsv} text "{{0.label_one_letter_code}}{{0.number}}" height 1.5
color /C,F #d9b38c
nucleotides /C,F fill
show ligand
style ligand stick
color ligand byelement
lighting soft
graphics silhouettes true
view
# >>> rotate, then: save Q9V2G1_IG.png width 2400 height 1800 supersample 3
""")
(H/"view.pml").write_text(f"""# Q9V2G1 IG on tRNA-bound 5WT1.  Open: pymol view.pml
load Q9V2G1_IG.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set cartoon_ring_mode,3
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select igtop, chain A and resi {"+".join(map(str,topres))}
color 0xcb181d, igtop
show sticks, igtop
util.cnc igtop
label igtop and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,14
select rna, chain C+F
show cartoon, rna
color wheat, rna
select cof, resn SAH
show sticks, cof
color grey60, cof
util.cnc cof
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png Q9V2G1_IG.png, dpi=300
""")
print("top-10 IG residues:", ", ".join(f"{r.residue}{int(r.position)}" for _,r in top.iterrows()))
print("wrote Q9V2G1_IG.pdb + view.pml + view.cxc")
