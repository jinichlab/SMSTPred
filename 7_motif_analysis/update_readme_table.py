#!/usr/bin/env python3
"""Regenerate the README structure table from the evidence-based rebuild.

Reads structures/motif_placements.csv (best in-class carrier per class-specific
motif, with SAM+substrate confirmed bound and placement measured) and rewrites the
"mapped structure and structural placement" section of README.md: per-class tables
of Consensus | Structure (protein, experiment, substrate + cofactor) | Placement
(near substrate / near SAM / near substrate+SAM / distal, with distances).
"""
import re
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
PL = pd.read_csv(HERE / "structures" / "motif_placements.csv")
README = HERE / "README.md"
CLASSES = ["DNA", "RNA", "Protein", "Metabolite"]

# concise enzyme labels keyed by UniProt accession (PDB disambiguates further)
SHORT = {
    "P05102": "M.HhaI", "Q9AXT8": "CMT3 chromomethylase", "P14385": "M.TaqI",
    "P08957": "EcoKI (Type I R-M)", "Q5C9L7": "norcoclaurine 6-OMT",
    "Q0CCY5": "GedA OMT", "A0A0B5L781": "MpaG' OMT", "I3V6A7": "scoulerine 9-OMT",
    "Q9H867": "METTL21D/VCPKMT", "Q63009": "PRMT1", "P0ACC1": "PrmC/HemK",
    "Q84BQ9": "PrmA (L11 MTase)", "Q8TEA1": "NSUN6", "Q58293": "aTrm5",
    "Q9UET6": "TrmJ/FTSJ1",
}


def short(row):
    return SHORT.get(row.accession, str(row.protein).split(" (")[0][:42])


def substrate(row):
    return {"DNA": "DNA", "RNA": "RNA", "Protein": "protein substrate"}.get(
        row["class"], f"small-molecule substrate (`{row.substrate}`)")


def place_cell(row):
    p = row.placement
    return f"**{p}** — SAM {row.d_SAM}, sub {row.d_substrate} Å"


def build():
    L = ["### The class-specific motifs — mapped structure and structural placement", ""]
    L += [
     "Of the 32 discovered motifs, **26 are class-specific** (own-class carrier rate ≥ 2× the",
     "pooled-background rate in the FIMO cross-scan over all 118,816 proteins; the other 6 are",
     "pan-methyltransferase SAM-binding / Rossmann-core motifs — see §2.1). For each of the 26 we",
     "select, **from FIMO evidence rather than by hand**, the best in-class carrier whose",
     "experimental structure has **both** a SAM/SAH cofactor and the substrate bound, locate the",
     "motif by its FIMO `matched_sequence` (author numbering), and measure the minimum",
     "heavy-atom distance from the motif to the cofactor and to the substrate. A motif is **near**",
     "a ligand at **≤ 4 Å** (direct van der Waals / H-bond contact); placement is **near substrate**,",
     "**near SAM**, **near substrate/SAM** (active site), or **distal** (> 4 Å from *both*).",
     "A geometry gate rejects structures where the ligand sits on a different",
     "chain/subunit (e.g. a ribosome-assembly complex), falling back to the next carrier.", "",
     "**#** is a sequential index within class; **Motif** is the original STREME ID (gaps are the",
     "6 non-specific motifs removed above). Per-motif PDB, chain, range and distances are in",
     "`motif_placements.csv`; the folder for each is `structures/<Class>/<class>_motif<STREME#>_<PDB>/`.",
     "See `structures/rebuild_from_carriers.py` for the pipeline.", "",
    ]
    for cls in CLASSES:
        g = PL[PL["class"] == cls].copy()
        g["num"] = g.motif.str.split("-").str[1].astype(int)
        g = g.sort_values("num")
        if g.empty:
            continue
        L.append(f"**{cls}**")
        L.append("")
        L.append("| # | Motif | Consensus | Structure | Placement |")
        L.append("|---|---|---|---|---|")
        for i, (_, r) in enumerate(g.iterrows(), 1):
            struct = f"{short(r)} `{r.pdb}` — {r.method}; **{r.cofactor} + {substrate(r)}**"
            L.append(f"| {i} | {r.motif} | `{r.consensus}` | {struct} | {place_cell(r)} |")
        L.append("")

    npl = PL.placement.value_counts()
    distal = PL[PL.placement == "distal"].sort_values("d_substrate")
    dnames = ", ".join(f"{r.motif} (`{r.consensus}`)" for _, r in distal.iterrows())
    L += [
     "**Summary.** All **26** class-specific motifs map onto an **experimental** structure of a",
     f"real MTase that carries the motif (X-ray; cryo-EM for the FTSJ1 panel), with **both** SAM/SAH",
     "and substrate bound and the placement measured directly — no homology models, no AlphaFold-only",
     "mappings, no hand-picked ranges. The motifs are overwhelmingly **active-site elements**:",
     f"{int(npl.get('near substrate/SAM',0))} sit at the substrate–cofactor junction (both ligands ≤ 4 Å), "
     f"{int(npl.get('near substrate',0))} contact only the substrate, {int(npl.get('near SAM',0))} only the SAM pocket,",
     f"and **{int(npl.get('distal',0))}** are **distal** — > 4 Å from *both* ligands: {dnames}.",
     "Of these, Metabolite-5, RNA-7 and RNA-6 sit well off the cofactor (> 9 Å from SAM) on the",
     "substrate-recognition face, while DNA-4 lies ~5–6 Å from both in the active-site cleft",
     "(junction-adjacent, in second shell rather than direct contact). The 6 non-specific",
     "motifs (RNA-3, Protein-1, Metabolite-1/2/3/4) are excluded as pan-MTase SAM-binding / Rossmann-core",
     "elements shared across classes.", "",
    ]
    return "\n".join(L)


def main():
    text = README.read_text()
    new = build()
    # replace the whole structure-table section up to (but not including) "## 4."
    pat = re.compile(r"### The (?:32 motifs|class-specific motifs).*?(?=\n---\n\n## 4\. Motif prevalence)", re.S)
    if not pat.search(text):
        raise SystemExit("section markers not found — check README")
    text = pat.sub(new + "\n", text)
    README.write_text(text)
    print(f"[updated] {README}  (+{len(new.splitlines())} lines)")
    print(new[:1500])


if __name__ == "__main__":
    main()
