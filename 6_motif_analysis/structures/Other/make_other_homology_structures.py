#!/usr/bin/env python3
"""
HOMOLOGY-mapped per-motif structures for the Other (small-molecule) motifs.

Our matched Other-MTases have no SAM-bound structures, so each motif's STREME PWM
was FIMO-matched onto a classic SAH-bound plant O-MTase co-crystal (the motif is a
conserved feature present in the homolog). Positions are the homolog's own author
numbering (located by the matched substring).

  1KYZ  caffeic-acid 3-O-MTase + SAH + FER (ferulic acid)  -> motifs 1,5,6,7,8
  1M6E  salicylic-acid carboxyl MTase + SAH + SAL          -> motif 4
  3C3Y  PFOMT phenylpropanoid/flavonoid O-MTase + SAH       -> motif 2 (no substrate)

Each PDB keeps the mapped protein chain + SAH cofactor + substrate ligand (waters/
ions/other chains stripped). view scripts: grey protein, motif coloured + sticks
(CPK), SAH grey sticks, substrate gold sticks, H-bond dashes.  *** HOMOLOG, not the
original matched protein -- label figures accordingly. ***
"""
import csv
from pathlib import Path
O = Path(__file__).resolve().parent
HOM = O / "homology"

# n, pdb, chain, substrate_resn, a_s, a_e, hex, folder, role
MOTIFS = [
    (1, "1KYZ", "E", "FER", 285, 295, "2171b5", "other_motif1", "caffeic-acid OMT homolog (SAH+ferulate)"),
    (3, "1KYZ", "E", "FER", 204, 215, "ce1256", "other_motif3_SAMloop", "caffeic-acid OMT homolog (SAH+ferulate) SAM Motif I loop"),
    (5, "1KYZ", "E", "FER", 326, 338, "e6550d", "other_motif5", "caffeic-acid OMT homolog (SAH+ferulate)"),
    (6, "1KYZ", "E", "FER", 268, 280, "238b45", "other_motif6", "caffeic-acid OMT homolog (SAH+ferulate)"),
    (7, "1KYZ", "E", "FER", 250, 265, "00868b", "other_motif7", "caffeic-acid OMT homolog (SAH+ferulate)"),
    (8, "1KYZ", "E", "FER", 227, 237, "8c564b", "other_motif8", "caffeic-acid OMT homolog (SAH+ferulate)"),
    (4, "1M6E", "X", "SAL", 128, 145, "cb181d", "other_motif4", "salicylate MTase homolog (SAH+salicylate)"),
    (2, "3C3Y", "A", None,  99, 111,  "6a51a3", "other_motif2", "PFOMT homolog (SAH, no substrate)"),
]


def subset(pdb, chain, substrate):
    keep_het = {"SAH"} | ({substrate} if substrate else set())
    out = []
    for ln in open(HOM / f"{pdb}.pdb"):
        rec = ln[:6]
        if rec == "CRYST1": out.append(ln); continue
        if rec not in ("ATOM  ", "HETATM", "TER   "): continue
        ch = ln[21]; rn3 = ln[17:20].strip()
        if rec.strip() in ("ATOM", "TER") and ch == chain:
            out.append(ln)
        elif rec == "HETATM" and rn3 in keep_het:
            out.append(ln)
    out.append("END\n")
    return out


def scripts(folder, chain, substrate, a_s, a_e, hx, role, n):
    sub_cx = f"color :{substrate} gold\n" if substrate else ""
    sub_pml = f"color yellow, resn {substrate}\nutil.cnc resn {substrate}\n" if substrate else ""
    restrict = "ligand"
    cxc = f"""# Other {role} [HOMOLOG].  Open: chimerax view.cxc
open {folder}.pdb
set bgColor white
hide atoms
cartoon
color /{chain} #cccccc
color /{chain}:{a_s}-{a_e} #{hx}
show /{chain}:{a_s}-{a_e} atoms
style /{chain}:{a_s}-{a_e} stick
color /{chain}:{a_s}-{a_e} byhetero
label /{chain}:{a_s}-{a_e} text "{{0.label_one_letter_code}}{{0.number}}" height 1.2
show ligand
style ligand stick
color ligand byelement
{sub_cx}hbonds /{chain}:{a_s}-{a_e} restrict ligand reveal true color #303030 dashes 8 radius 0.05
2dlabel text "Other motif{n} (homolog): {role}" xpos .03 ypos .04 size 13
lighting soft
graphics silhouettes true
view
# >>> rotate, then: save {folder}.png width 2400 height 1800 supersample 3
"""
    pml = f"""# Other {role} [HOMOLOG].  Open: pymol view.pml
load {folder}.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain {chain}
color grey80, chain {chain}
select motif, chain {chain} and resi {a_s}-{a_e}
color 0x{hx}, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
show sticks, resn SAH
color grey50, resn SAH
util.cnc resn SAH
{sub_pml}distance contacts, motif, (resn SAH{' or resn '+substrate if substrate else ''}), 3.6, mode=2
hide labels, contacts
color black, contacts
orient chain {chain}
deselect
# >>> rotate, then: ray 2400,1800; png {folder}.png, dpi=300
"""
    return pml, cxc


rows = []
for n, pdb, ch, sub, a_s, a_e, hx, folder, role in MOTIFS:
    d = O / folder; d.mkdir(exist_ok=True)
    lines = subset(pdb, ch, sub)
    (d / f"{folder}.pdb").write_text("".join(lines))
    pml, cxc = scripts(folder, ch, sub, a_s, a_e, hx, role, n)
    (d / "view.pml").write_text(pml); (d / "view.cxc").write_text(cxc)
    natom = sum(1 for x in lines if x[:4] == "ATOM")
    print(f"  motif {n} [{pdb} /{ch}:{a_s}-{a_e} + SAH{('+'+sub) if sub else ''}] -> {folder}/ ({natom} atoms)")
    rows.append({"motif": n, "folder": folder, "structure": pdb, "type": "homology",
                 "motif_residues": f"/{ch}:{a_s}-{a_e}", "cofactor": "SAH",
                 "substrate": sub or "-", "role": role})
with open(O / "other_homology_structures.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("[done]")
