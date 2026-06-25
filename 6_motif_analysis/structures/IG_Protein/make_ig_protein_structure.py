#!/usr/bin/env python3
"""IG of Q54S83 (NDUFAF7, no-motif protein-arginine MTase) mapped onto its apo
structure 5ZTZ. IG -> B-factor (chain A). Keep chain A only (drop extra copies B/C,
waters). view: grey protein, top-10 IG residues red sticks + labels (apo: no
substrate/cofactor)."""
from pathlib import Path
import pandas as pd
H=Path(__file__).resolve().parent
ig=pd.read_csv(H.parent.parent/"attribution"/"ig_q54s83"/"Q54S83__Protein_residue_attr.csv")
igmap=dict(zip(ig.position.astype(int),ig.ig.astype(float)))
top=ig.sort_values("ig",ascending=False).head(10); topres=list(top.position.astype(int))
out=[]
for ln in open(H/"5ZTZ.pdb"):
    rec=ln[:6]
    if rec=="CRYST1": out.append(ln); continue
    if rec in ("ATOM  ","TER   ") and ln[21]=="A":
        if rec=="ATOM  ":
            try: r=int(ln[22:26]); ln=ln[:60]+f"{igmap.get(r,0.0)*100:6.2f}"+ln[66:]
            except: pass
        out.append(ln)
out.append("END\n")
(H/"Q54S83_IG.pdb").write_text("".join(out))
sel=",".join(map(str,topres)); selp="+".join(map(str,topres))
(H/"view.cxc").write_text(f"""# Q54S83 (NDUFAF7, no-motif protein-MTase) IG on apo 5ZTZ.  Open: chimerax view.cxc
open Q54S83_IG.pdb
set bgColor white
hide atoms
cartoon
color /A #cccccc
color /A:{sel} #cb181d
show /A:{sel} atoms
style /A:{sel} stick
color /A:{sel} byhetero
label /A:{sel} text "{{0.label_one_letter_code}}{{0.number}}" height 1.5
lighting soft
graphics silhouettes true
view
# >>> rotate, then: save Q54S83_IG.png width 2400 height 1800 supersample 3
""")
(H/"view.pml").write_text(f"""# Q54S83 IG on apo 5ZTZ.  Open: pymol view.pml
load Q54S83_IG.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select igtop, chain A and resi {selp}
color 0xcb181d, igtop
show sticks, igtop
util.cnc igtop
label igtop and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,14
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png Q54S83_IG.png, dpi=300
""")
print("top-10 IG:", ", ".join(f"{r.residue}{int(r.position)}" for _,r in top.iterrows()))
print("wrote Q54S83_IG.pdb + view scripts")
