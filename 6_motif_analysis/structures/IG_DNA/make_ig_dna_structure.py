#!/usr/bin/env python3
"""
Map the per-residue Integrated-Gradients importance of T4 Dam (UniProt P04392 --
the DNA-adenine MTase we ran IG on) onto its DNA-bound structure 1YF3.

The IG value (signed, logit-space; + drives the DNA call) is written into the
B-factor column of protein chain A, so the protein can be coloured as an IG
heat-map (blue=away, white=0, red=toward-DNA-call). Output keeps chain A + the DNA
duplex (C,D) + the AdoHcy cofactor (SAH); waters/glycerol/2nd protomer removed.
1YF3 author numbering == UniProt numbering (verified), so IG positions map directly.

Output (structures/IG_DNA/):
  P04392_IG.pdb         chain A (B-factor=IG*100) + DNA + cofactor
  view.pml / view.cxc   colour-by-IG + DNA + cofactor + top-IG residues as sticks
"""
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
IG_CSV = HERE.parent.parent / "attribution" / "ig_output_logit" / "P04392__DNA_residue_attr.csv"
PDB = HERE / "1YF3.pdb"
KEEP_CHAINS = {"A", "C", "D"}          # protein A + DNA C,D
DROP_HET = {"HOH", "GOL"}


def main():
    ig = pd.read_csv(IG_CSV)
    ig_by_res = dict(zip(ig["position"].astype(int), ig["ig"].astype(float)))
    vmax = max(abs(v) for v in ig_by_res.values()) * 100
    top = ig.sort_values("ig", ascending=False).head(10)
    top_res = list(top["position"].astype(int))
    top_lbl = [f"{r.residue}{int(r.position)}" for _, r in top.iterrows()]

    out = []
    for ln in open(PDB):
        rec = ln[:6]
        if rec == "CRYST1":
            out.append(ln); continue
        if rec in ("ATOM  ", "HETATM", "TER   "):
            ch = ln[21]; resname = ln[17:20].strip()
            if ch not in KEEP_CHAINS:
                continue
            if rec == "HETATM" and resname in DROP_HET:
                continue
            if rec == "ATOM  " and ch == "A":     # write IG into B-factor
                try:
                    rn = int(ln[22:26])
                except ValueError:
                    out.append(ln); continue
                b = ig_by_res.get(rn, 0.0) * 100
                ln = ln[:60] + f"{b:6.2f}" + ln[66:]
            elif rec in ("ATOM  ", "HETATM"):       # DNA / cofactor: neutral B
                ln = ln[:60] + f"{0.0:6.2f}" + ln[66:]
            out.append(ln)
    out.append("END\n")
    (HERE / "P04392_IG.pdb").write_text("".join(out))
    print(f"wrote P04392_IG.pdb  | IG range +/-{vmax:.1f} (B-factor*100)")
    print("top-10 IG residues:", ", ".join(top_lbl))

    sel = "+".join(str(r) for r in top_res)
    # ---- PyMOL ----
    (HERE / "view.pml").write_text(f"""# T4 Dam (P04392) IG on DNA-bound 1YF3.  Open: pymol view.pml
load P04392_IG.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set ray_shadows, 0
# whole protein one colour (grey)
show cartoon, chain A
color grey80, chain A
# top-IG residues highlighted as red sticks (CPK heteroatoms) + labels
select igtop, chain A and resi {sel}
color 0xcb181d, igtop
show sticks, igtop
util.cnc igtop
label igtop and name CA, "%s%s" % (one_letter[resn], resi)
set label_size, 14
# DNA duplex (tan slabs)
show cartoon, chain C+D
color wheat, chain C+D
# cofactor AdoHcy
select cof, resn SAH
show sticks, cof
color grey50, cof
util.cnc cof
orient chain A
deselect
# >>> rotate, then:  ray 2400,1800;  png P04392_IG.png, dpi=300
""")
    # ---- ChimeraX ----
    cxc_sel = ",".join(str(r) for r in top_res)
    (HERE / "view.cxc").write_text(f"""# T4 Dam (P04392) IG on DNA-bound 1YF3.  Open: chimerax view.cxc
open P04392_IG.pdb
set bgColor white
hide atoms
cartoon
color /A #cccccc
color /A:{cxc_sel} #cb181d
show /A:{cxc_sel} atoms
style /A:{cxc_sel} stick
color /A:{cxc_sel} byhetero
label /A:{cxc_sel} text "{{0.label_one_letter_code}}{{0.number}}" height 1.5
color /C,D #d9b38c
nucleotides /C,D fill
show ligand
style ligand stick
color ligand byelement
lighting soft
graphics silhouettes true
view
# >>> rotate, then:  save P04392_IG.png width 2400 height 1800 supersample 3
""")
    print("wrote view.pml + view.cxc")


if __name__ == "__main__":
    main()
