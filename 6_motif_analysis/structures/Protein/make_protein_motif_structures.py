#!/usr/bin/env python3
"""
Per-motif PDB + pre-styled view scripts for the Protein STREME motifs.

  2NXN  PrmA (L11 MTase) chain A + L11 protein SUBSTRATE (chain B)  -> motifs 1, 7
  4GQB  PRMT5 chain A + MEP50 partner (B) + peptide (C) + cofactor   -> motifs 4, 6
  AlphaFold (v6) monomer models                                     -> motifs 2,3,5,8
        m2 Q6C4Y5 (SET lysine MTase), m3 A7GSD9 (O-MTase),
        m5 Q8DHV7 (HemK/PrmC), m8 Q9A9T7 (HemK, NPPY catalytic)

Complexes keep the whole MTase chain + partner protein(s) + cofactor (waters/ions
stripped); AF models keep the monomer. view.pml/view.cxc load pre-styled: MTase grey,
motif coloured + sticks (CPK heteroatoms) + label, partner protein wheat cartoon,
cofactor sticks, H-bond dashes -- open, rotate, save.
"""
import csv
from pathlib import Path
P = Path(__file__).resolve().parent

# n, pdb, mtase_chain, partner_chains, cofactor, up_s, up_e, hex, folder, role
MOTIFS = [
    (1, "2NXN", "A", ["B"], None, 126, 141, "2171b5", "protein_motif1_PrmA_SAMloop", "PrmA SAM Motif I loop (+L11 substrate)"),
    (7, "2NXN", "A", ["B"], None, 98, 114,  "00868b", "protein_motif7_PrmA",          "PrmA region (+L11 substrate)"),
    (4, "4GQB", "A", ["B","C"], "0XU", 429, 445, "cb181d", "protein_motif4_PRMT5",     "PRMT5 region (+MEP50 partner)"),
    (6, "4GQB", "A", ["B","C"], "0XU", 45, 53,   "ce1256", "protein_motif6_PRMT5",     "PRMT5 region (+MEP50 partner)"),
    (2, "AF-Q6C4Y5", "A", [], None, 336, 350, "6a51a3", "protein_motif2_SET",   "SET lysine-MTase region (AlphaFold)"),
    (3, "AF-A7GSD9", "A", [], None, 196, 203, "238b45", "protein_motif3",       "O-MTase region (AlphaFold)"),
    (5, "AF-Q8DHV7", "A", [], None, 10, 19,   "e6550d", "protein_motif5_HemK",  "HemK/PrmC region (AlphaFold)"),
    (8, "AF-Q9A9T7", "A", [], None, 188, 195, "8c564b", "protein_motif8_HemK_NPPY", "HemK catalytic NPPY (AlphaFold)"),
]
DROPHET = {"HOH","GOL","EDO","ACT","PO4","SO4","NA","CL","K","MG","ZN","ACE"}


def subset(pdb, mt, partners, cof):
    keep = set([mt]) | set(partners)
    out = []
    for ln in open(P / f"{pdb}.pdb"):
        rec = ln[:6]
        if rec == "CRYST1": out.append(ln); continue
        if rec not in ("ATOM  ", "HETATM", "TER   "): continue
        ch = ln[21]; rn3 = ln[17:20].strip()
        if rec.strip() in ("ATOM", "TER") and ch in keep:
            out.append(ln)
        elif rec == "HETATM" and cof and rn3 == cof:
            out.append(ln)
    out.append("END\n")
    return out


def scripts(folder, mt, partners, cof, a_s, a_e, hx, role, n):
    has_partner = bool(partners)
    pml_partner = (f"select partner, chain {'+'.join(partners)}\nshow cartoon, partner\ncolor wheat, partner\n"
                   if has_partner else "")
    pml_cof = ("select cof, resn %s\nshow sticks, cof\ncolor grey50, cof\nutil.cnc cof\n" % cof) if cof else ""
    pml_hb = (f"distance contacts, motif, (chain {'+'.join(partners)}{' or resn '+cof if cof else ''}), 3.6, mode=2\nhide labels, contacts\ncolor black, contacts\n"
              if (has_partner or cof) else "")
    pml = f"""# Protein {role}.  Open: pymol view.pml
load {folder}.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain {mt}
color grey80, chain {mt}
select motif, chain {mt} and resi {a_s}-{a_e}
color 0x{hx}, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
{pml_partner}{pml_cof}{pml_hb}orient chain {mt}
deselect
# >>> rotate, then: ray 2400,1800; png {folder}.png, dpi=300
"""
    pc_cx = ",".join(f"/{c}" for c in partners)
    cxc_partner = (f"color {pc_cx} #d9b38c\n") if has_partner else ""
    cxc_cof = "show ligand\nstyle ligand stick\ncolor ligand byelement\n" if cof else ""
    restrict = " | ".join(([pc_cx] if has_partner else []) + (["ligand"] if cof else []))
    cxc_hb = (f"hbonds /{mt}:{a_s}-{a_e} restrict {restrict} reveal true color #303030 dashes 8 radius 0.05\n"
              if restrict else "")
    cxc = f"""# Protein {role}.  Open: chimerax view.cxc
open {folder}.pdb
set bgColor white
hide atoms
cartoon
color /{mt} #cccccc
color /{mt}:{a_s}-{a_e} #{hx}
show /{mt}:{a_s}-{a_e} atoms
style /{mt}:{a_s}-{a_e} stick
color /{mt}:{a_s}-{a_e} byhetero
label /{mt}:{a_s}-{a_e} text "{{0.label_one_letter_code}}{{0.number}}" height 1.2
{cxc_partner}{cxc_cof}{cxc_hb}2dlabel text "Protein motif{n}: {role}" xpos .03 ypos .04 size 14
lighting soft
graphics silhouettes true
view
# >>> rotate, then: save {folder}.png width 2400 height 1800 supersample 3
"""
    return pml, cxc


rows = []
for n, pdb, mt, partners, cof, us, ue, hx, folder, role in MOTIFS:
    d = P / folder; d.mkdir(exist_ok=True)
    lines = subset(pdb, mt, partners, cof)
    natom = sum(1 for x in lines if x[:4] == "ATOM")
    (d / f"{folder}.pdb").write_text("".join(lines))
    pml, cxc = scripts(folder, mt, partners, cof, us, ue, hx, role, n)
    (d / "view.pml").write_text(pml); (d / "view.cxc").write_text(cxc)
    typ = "complex" if not pdb.startswith("AF") else "AlphaFold"
    print(f"  motif {n} [{pdb} /{mt}:{us}-{ue}] {typ} -> {folder}/  ({natom} atoms)")
    rows.append({"motif": n, "folder": folder, "structure": pdb, "type": typ,
                 "motif_residues": f"/{mt}:{us}-{ue}", "partner": "+".join(partners) or "-",
                 "cofactor": cof or "-", "role": role})
with open(P / "protein_motif_structures.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("[done]")
