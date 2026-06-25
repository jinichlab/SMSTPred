#!/usr/bin/env python3
"""
Write one PDB per DNA STREME motif for manual figure-making.

Each output keeps: the motif protein region (+/-8 flanking residues, chain A) +
the FULL DNA duplex (chains B,C) + the cofactor (SAH for HhaI / NEA AdoMet-analog
for TaqI) + the flipped/target base. The bulk protein (and TaqI's 2nd complex,
chains D/E/F) and waters are removed. Coordinates/numbering are taken verbatim
from the deposited structure, so the exact motif residues colour cleanly
(e.g. `/A:76-88`). A manifest lists each motif's exact residues, role, and contacts.

Output: structures/DNA/<folder>/<folder>.pdb  +  dna_motif_structures.csv
"""
import csv
from pathlib import Path

DNA = Path(__file__).resolve().parent
FLANK = 8
DNA_CHAINS = {"B", "C"}

# n, pdb, s, e, folder, consensus, structure-derived role, measured contacts
MOTIFS = [
    (1, "1MHT", 76, 88,  "motif1_PCQ",      "IGGPPCQGFSRAG",      "Motif IV catalytic PCQ (Cys81)",     "DNA flipped-base 1.8A + cofactor 3.6A"),
    (2, "1MHT", 17, 24,  "motif2_SAMloop",  "PFAGSGGT",           "Motif I SAM glycine loop",           "cofactor 2.9A (DNA 12.7A)"),
    (3, "1MHT", 157, 167,"motif3_TRDbasic", "YGVPQRRERVF",        "TRD DNA-contacting basic loop",      "DNA 2.8A / base 2.6A"),
    (4, "1MHT", 274, 291,"motif4_TRD",      "LTPREAARLQGFPDDYTF", "cofactor-proximal (NOT DNA/TRD)",    "cofactor 3.8A (DNA 9.6A far)"),
    (5, "1MHT", 300, 315,"motif5_TRD",      "KQIGNAVPPPVARA",     "DNA-cofactor junction",              "DNA 5.2 / base 3.4 / cofactor 3.2A"),
    (6, "1MHT", 116, 124,"motif6_GluVI",    "FLLENVPGL",          "Motif VI catalytic Glu119",          "DNA 3.8A / base 2.8A"),
    (7, "1G38", 103, 110,"motif7_NPPY",     "YTBPPYNT",           "Motif IV catalytic NPPY (Asn105)",   "DNA 2.9A + cofactor 2.9A"),
    (8, "1G38", 19, 27,  "motif8_SAMregion","GZFFTPREI",          "DNA-cofactor junction (NOT SAM-only)","DNA 3.3A + cofactor 3.8A"),
]


def subset(pdb, s_lo, s_hi):
    lines = []
    for ln in open(DNA / f"{pdb}.pdb"):
        rec = ln[:6]
        if rec == "CRYST1":
            lines.append(ln); continue
        if rec in ("ATOM  ", "HETATM", "TER   "):
            chain = ln[21]
            if pdb == "1G38" and chain in {"D", "E", "F"}:
                continue
            resname = ln[17:20].strip()
            if rec.startswith("ATOM") or rec == "TER   ":
                if chain == "A":
                    try:
                        rs = int(ln[22:26])
                    except ValueError:
                        continue
                    if s_lo <= rs <= s_hi:
                        lines.append(ln)
                elif chain in DNA_CHAINS:
                    lines.append(ln)
            else:  # HETATM: keep cofactor / modified bases, drop water
                if resname != "HOH" and chain in {"A", "B", "C"}:
                    lines.append(ln)
    lines.append("END\n")
    return lines


def main():
    rows = []
    for n, pdb, s, e, folder, cons, role, contacts in MOTIFS:
        d = DNA / folder
        d.mkdir(exist_ok=True)
        out = d / f"{folder}.pdb"
        lines = subset(pdb, s - FLANK, e + FLANK)
        out.write_text("".join(lines))
        n_atom = sum(1 for x in lines if x[:4] == "ATOM")
        n_het = sum(1 for x in lines if x[:6] == "HETATM")
        cof = "SAH" if pdb == "1MHT" else "NEA"
        print(f"  motif {n} [{pdb} /A:{s}-{e}] -> {folder}/{folder}.pdb  "
              f"({n_atom} atom, {n_het} het)")
        rows.append({"motif": n, "folder": folder, "pdb": out.name,
                     "source_structure": pdb, "motif_residues": f"/A:{s}-{e}",
                     "consensus": cons, "cofactor": cof, "role": role,
                     "contacts": contacts})
    with open(DNA / "dna_motif_structures.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"[done] 8 motif PDBs + dna_motif_structures.csv")


if __name__ == "__main__":
    main()
