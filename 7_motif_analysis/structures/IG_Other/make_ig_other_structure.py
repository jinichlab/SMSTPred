#!/usr/bin/env python3
"""IG of P0DPA9 (PsiM psilocybin synthase, no-motif Other MTase) on 8PB3
(SAH + norbaeocystin substrate XP6). IG->B-factor (chain A); keep chain A + SAH +
XP6. view: grey protein, top-10 IG red sticks+labels, SAH grey sticks, XP6 gold."""
from pathlib import Path
import pandas as pd
H=Path(__file__).resolve().parent
ig=pd.read_csv(H.parent.parent/"attribution"/"ig_p0dpa9"/"P0DPA9__Other_residue_attr.csv")
igmap=dict(zip(ig.position.astype(int),ig.ig.astype(float)))
top=ig.sort_values("ig",ascending=False).head(10); topres=list(top.position.astype(int))
KEEPHET={"SAH","XP6"}
out=[]
for ln in open(H/"8PB3.pdb"):
    rec=ln[:6]
    if rec=="CRYST1": out.append(ln); continue
    if rec in ("ATOM  ","TER   ") and ln[21]=="A":
        if rec=="ATOM  ":
            try: r=int(ln[22:26]); ln=ln[:60]+f"{igmap.get(r,0.0)*100:6.2f}"+ln[66:]
            except: pass
        out.append(ln)
    elif rec=="HETATM" and ln[21]=="A" and ln[17:20].strip() in KEEPHET:
        out.append(ln[:60]+f"{0.0:6.2f}"+ln[66:])
out.append("END\n")
(H/"P0DPA9_IG.pdb").write_text("".join(out))
sel=",".join(map(str,topres)); selp="+".join(map(str,topres))
(H/"view.cxc").write_text(f"""# P0DPA9 (PsiM, no-motif Other MTase) IG on 8PB3 (SAH+substrate).  Open: chimerax view.cxc
open P0DPA9_IG.pdb
set bgColor white
hide atoms
cartoon
color /A #cccccc
color /A:{sel} #cb181d
show /A:{sel} atoms
style /A:{sel} stick
color /A:{sel} byhetero
label /A:{sel} text "{{0.label_one_letter_code}}{{0.number}}" height 1.5
show ligand
style ligand stick
color ligand byelement
color :XP6 gold
lighting soft
graphics silhouettes true
view
# >>> rotate, then: save P0DPA9_IG.png width 2400 height 1800 supersample 3
""")
(H/"view.pml").write_text(f"""# P0DPA9 IG on 8PB3 (SAH+substrate).  Open: pymol view.pml
load P0DPA9_IG.pdb, m
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
show sticks, resn SAH+XP6
util.cnc resn SAH+XP6
color yellow, resn XP6 and elem C
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png P0DPA9_IG.png, dpi=300
""")
print("top-10 IG:", ", ".join(f"{r.residue}{int(r.position)}" for _,r in top.iterrows()))
