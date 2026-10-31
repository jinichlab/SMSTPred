#!/usr/bin/env python3
"""
Make every motif folder self-contained: write a FULL-structure PDB + matching
view.cxc / view.pml that show the WHOLE protein (grey ribbon, no side chains)
with the motif highlighted (colour + CPK sticks) and substrate/cofactor in place.

Someone can pull a single motif folder (the .pdb + a view file) and reproduce the
context figure.  Replaces the earlier mix of trimmed (DNA, RNA-complex) and
full-chain (Protein, Metabolite, RNA-apo) motif PDBs with one consistent output.

Per motif we keep: the full MTase chain + substrate chains (DNA/RNA/partner
protein) + substrate/cofactor ligands; waters, ions and unrelated chains dropped.
Config (source PDB, chain, substrate kind, residue range, colour) is lifted from
the per-class make_*/render_* scripts.  Structure (author) residue numbering.
"""
import shutil
from pathlib import Path

S = Path(__file__).resolve().parent
WATER_ION = {"HOH", "GOL", "EDO", "ACT", "PO4", "SO4", "NA", "CL", "K", "MG",
             "ZN", "CA", "MN", "NI", "FE", "BME", "DMS", "IOD", "BR", "ACE",
             "FMT", "ADP", "EPE", "HF2", "UNX"}
KEEP_LIGS = {"SAH", "SAM", "0XU", "NEA", "C36", "FER", "SAL", "MLZ", "SFG",
             "QVR", "PPY", "M3L"}

# The SAM cofactor and its analogs get distinct warm-brown carbons, so the cofactor
# never reads as "protein grey" and never collides with the gold substrate ligand.
# (ALY / MLZ / C36 / 6MA are modified residues, not cofactors -- left alone.)
COF_LIGS = ("0XU", "NEA", "QVR", "SAH", "SAM", "SFG")
COF_HEX = "#b5651d"                  # warm brown-orange
COF_RGB = "[0.710, 0.396, 0.114]"    # same colour, PyMOL float triplet

# cls, folder, src, src_subdir, chain, sub_chains, sub_kind, sub_lig, show_cof, rng, hex, caption
M = [
 ("DNA","motif1_PCQ","1MHT","","A",["B","C"],"nucleic",None,True,"76-88","cb181d","motif1 IGGPPCQGFSRAG - catalytic PCQ (Cys81)"),
 ("DNA","motif2_SAMloop","1MHT","","A",["B","C"],"nucleic",None,True,"17-24","2171b5","motif2 - Motif I SAM glycine loop"),
 ("DNA","motif3_TRDbasic","1MHT","","A",["B","C"],"nucleic",None,True,"157-167","238b45","motif3 - TRD DNA-contacting basic loop"),
 ("DNA","motif4_TRD","1MHT","","A",["B","C"],"nucleic",None,True,"274-291","6a51a3","motif4 - target-recognition domain"),
 ("DNA","motif5_TRD","1MHT","","A",["B","C"],"nucleic",None,True,"300-315","ce1256","motif5 - target-recognition domain"),
 ("DNA","motif6_GluVI","1MHT","","A",["B","C"],"nucleic",None,True,"116-124","e6550d","motif6 - Motif VI catalytic Glu119"),
 ("DNA","motif7_NPPY","1G38","","A",["B","C"],"nucleic",None,True,"103-110","cb181d","motif7 - Motif IV catalytic NPPY (Asn105)"),
 ("DNA","motif8_SAMregion","1G38","","A",["B","C"],"nucleic",None,True,"19-27","2171b5","motif8 - Motif I SAM region"),
 ("Protein","protein_motif1_PrmA_SAMloop","2NXN","","A",["B"],"partner",None,False,"126-141","2171b5","motif1 - PrmA SAM Motif I loop (+L11 substrate)"),
 ("Protein","protein_motif7_PrmA","2NXN","","A",["B"],"partner",None,False,"98-114","00868b","motif7 - PrmA region (+L11 substrate)"),
 ("Protein","protein_motif4_PRMT5","4GQB","","A",["B","C"],"partner",None,True,"429-445","cb181d","motif4 - PRMT5 region (+MEP50)"),
 ("Protein","protein_motif6_PRMT5","4GQB","","A",["B","C"],"partner",None,True,"45-53","ce1256","motif6 - PRMT5 region (+MEP50)"),
 # NOTE: Protein motif2/3/8 and Metabolite motif2 moved to make_context_views_upgrades.py --
 # each now uses a substrate-BOUND experimental structure found by scanning all FIMO
 # carriers (8HL7 METTL21D+p97, 2B3T PrmC+RF1, 4KIC + phenylpyruvate). Do not re-add.
 ("RNA","rna_motif2_METTL16","9U48","","A",["B"],"nucleic",None,True,"135-142","cb181d","motif2 - METTL16 catalytic region"),
 ("RNA","rna_motif3_SAMloop","2ZZM","","A",["B"],"nucleic",None,True,"198-212","2171b5","motif3 - Trm5 SAM DxGxGxG loop"),
 ("RNA","rna_motif6","2ZZM","","A",["B"],"nucleic",None,True,"45-55","6a51a3","motif6 - Trm5 region"),
 ("RNA","rna_motif7_basic","2ZZM","","A",["B"],"nucleic",None,True,"233-240","e6550d","motif7 - Trm5 basic RNA-binding patch"),
 # NOTE: RNA motif1/4/5/8 and Protein motif5 are built by make_context_views_upgrades.py
 # from substrate/cofactor-BOUND structures (NSUN6 5WWS, FtsJ 1EJ0, CARM1 7OKP) that the
 # cluster-member search found. Do not re-add them here (would overwrite with apo/AF).
 ("Metabolite","metabolite_motif1","1KYZ","homology","E",[],"ligand","FER",True,"285-295","2171b5","motif1 - caffeic-acid OMT homolog (SAH+ferulate)"),
 ("Metabolite","metabolite_motif3_SAMloop","1KYZ","homology","E",[],"ligand","FER",True,"204-215","ce1256","motif3 - SAM Motif I loop (homolog)"),
 ("Metabolite","metabolite_motif5","1KYZ","homology","E",[],"ligand","FER",True,"326-338","e6550d","motif5 - caffeic-acid OMT homolog"),
 ("Metabolite","metabolite_motif6","1KYZ","homology","E",[],"ligand","FER",True,"268-280","238b45","motif6 - caffeic-acid OMT homolog"),
 ("Metabolite","metabolite_motif7","1KYZ","homology","E",[],"ligand","FER",True,"250-265","00868b","motif7 - caffeic-acid OMT homolog"),
 ("Metabolite","metabolite_motif8","1KYZ","homology","E",[],"ligand","FER",True,"227-237","8c564b","motif8 - caffeic-acid OMT homolog"),
 ("Metabolite","metabolite_motif4","1M6E","homology","X",[],"ligand","SAL",True,"128-145","cb181d","motif4 - salicylate MTase homolog (SAH+salicylate)"),
]


def subset(src, keep_chains):
    out = []
    for ln in open(src):
        rec = ln[:6]
        if rec == "CRYST1":
            out.append(ln); continue
        if rec not in ("ATOM  ", "HETATM", "TER   "):
            continue
        ch = ln[21]; resn = ln[17:20].strip()
        if rec in ("ATOM  ", "TER   "):
            if ch in keep_chains:
                out.append(ln)
        else:  # HETATM: only keep cofactor/substrate bound to a kept chain
            if resn in WATER_ION:
                continue
            if ch in keep_chains:      # ligands are assigned to their protein's chain
                out.append(ln)
    out.append("END\n")
    return out


def cxc(folder, chain, sub_chains, kind, sub_lig, cof, rng, hx, cap):
    sub, tgt, hide_nuc = "", [], ""
    if kind == "nucleic":
        # colour the whole nucleic chains tan (standard + modified/HETATM nucleotides)
        sub = "".join(f"color /{c} #9c7a5b\n" for c in sub_chains)
        # hide the atom/stick display on the nucleic chains -> every nucleotide reads
        # as the uniform cartoon; modified/target bases no longer show as element sticks
        hide_nuc = "hide /" + ",".join(sub_chains) + " atoms\n"   # ChimeraX: /B,C
    elif kind == "partner":
        sub = "".join(f"color /{c} #d9b38c\n" for c in sub_chains); tgt += [f"/{c}" for c in sub_chains]
    cofblk = ""
    if cof:
        cofsel = ",".join(COF_LIGS)
        cofblk = ("show ligand\nstyle ligand stick\ncolor ligand byelement\n"
                  f"color :{cofsel} {COF_HEX}\ncolor :{cofsel} byhetero\n")
        if sub_lig:
            cofblk += f"color :{sub_lig} gold\ncolor :{sub_lig} byhetero\n"
        tgt.append("ligand")
    hb = ""
    if tgt:
        hb = (f"hbonds /{chain}:{rng} restrict {' | '.join(tgt)} reveal true "
              f"color #303030 dashes 8 radius 0.05\n")
    return f"""# {cap}.  Open:  chimerax view.cxc
open {folder}.pdb
set bgColor white
hide atoms
cartoon
color #1 #b0b0b0
{sub}{cofblk}{hide_nuc}color /{chain}:{rng} #{hx}
show /{chain}:{rng} atoms
style /{chain}:{rng} stick
color /{chain}:{rng} byhetero
{hb}2dlabel text "{cap}" xpos .03 ypos .04 size 14
lighting soft
graphics silhouettes true
view
# >>> rotate to taste, then:  save {folder}.png width 2400 height 1800 supersample 3
"""


def pml(folder, chain, sub_chains, kind, sub_lig, cof, rng, hx, cap):
    sub, tgt, cof_excl = "", [], ""
    if kind == "nucleic":
        # whole nucleic chains as uniform tan cartoon (standard + modified nucleotides);
        # no per-nucleotide sticks, so protein vs nucleic reads cleanly
        nsel = "+".join(sub_chains)
        sub = ("set_color dnatan, [0.612,0.478,0.357]\n"
               f"select nuc, chain {nsel}\nshow cartoon, nuc\nhide sticks, nuc\ncolor dnatan, nuc\n")
        cof_excl = f" and not chain {nsel}"      # keep modified/target bases out of the cofactor sticks
    elif kind == "partner":
        sel = "+".join(sub_chains)
        sub = f"select partner, chain {sel}\nshow cartoon, partner\ncolor wheat, partner\n"; tgt.append("partner")
    cofblk = ""
    if cof:
        resn = "+".join(COF_LIGS)
        cofblk = (f"select cof, (hetatm and not solvent{cof_excl})\n"
                  "show sticks, cof\ncolor grey60, cof\nutil.cnc cof\n"
                  f"set_color cofcarbon, {COF_RGB}\n"
                  f"select samcof, (cof and resn {resn})\n"
                  "color cofcarbon, samcof\nutil.cnc samcof\n")
        if sub_lig:
            cofblk += f"color yellow, resn {sub_lig}\nutil.cnc resn {sub_lig}\n"
        tgt.append("cof")
    hb = ""
    if tgt:
        hb = (f"distance contacts, motif, ({' or '.join(tgt)}), 3.6, mode=2\n"
              "hide labels, contacts\ncolor black, contacts\n")
    return f"""# {cap}.  Open:  pymol view.pml
load {folder}.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
set valence, 1
show cartoon, polymer.protein
color 0xb0b0b0, polymer.protein
select motif, chain {chain} and resi {rng}
color 0x{hx}, motif
show sticks, motif
util.cnc motif
{sub}{cofblk}{hb}set_title m, 1, "{cap}"
orient
deselect
# >>> rotate to taste, then:  ray 2400, 1800;  png {folder}.png, dpi=300
"""


def main():
    for cls, folder, src, subdir, chain, subch, kind, slig, cof, rng, hx, cap in M:
        d = S / cls / folder
        d.mkdir(parents=True, exist_ok=True)
        srcpath = S / cls / subdir / f"{src}.pdb"
        keep = {chain} | set(subch)
        (d / f"{folder}.pdb").write_text("".join(subset(srcpath, keep)))
        (d / "view.cxc").write_text(cxc(folder, chain, subch, kind, slig, cof, rng, hx, cap))
        (d / "view.pml").write_text(pml(folder, chain, subch, kind, slig, cof, rng, hx, cap))
        n_prot = sum(1 for _ in open(d / f"{folder}.pdb")
                     if _.startswith("ATOM") and _[21] == chain)
        print(f"  {cls}/{folder:30} {src} chain {chain} :{rng}  ({n_prot} MTase atoms)")


if __name__ == "__main__":
    main()
