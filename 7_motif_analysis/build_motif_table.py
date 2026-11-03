#!/usr/bin/env python3
"""Assemble the STREME motif table (Table 2): the class-specific motifs of every
substrate label (carrier rate in the motif's own class at least 2x the rate outside it),
with consensus, discovery stats, REAL per-protein coverage (fraction of that class's
members carrying it, from the FIMO-all run), fold enrichment and structural placement.

Inputs (this dir):
  streme_{cls}/streme.txt      8 STREME motifs per class (consensus, width, E-value)
  fimo_all_{cls}/fimo.tsv       FIMO hits of those motifs vs ALL class members
  class_{cls}.fasta             class member counts (denominator)
  python build_motif_table.py                 the archived table (this folder)
  python build_motif_table.py --dir <release>/motif_analysis    a refresh release: the STREME motifs of <dir>, coverage and fold from the FIMO cross-scan (motif_specificity_matrix.csv),
                              placement and distances from structure_search/motif_placements.csv (motif_structure_pick.py); writes motif_table.csv/.png into <dir>, with the PDB id of the structure
Outputs:
  motif_specificity_matrix.csv  fold enrichment (ratio_vs_bg) = own-class carrier rate / rate across all proteins outside the class
  structures/motif_placements.csv  structural placement of the 26 class-specific motifs (the other six are not mapped)
  motif_table.csv               one row per class-specific motif (26 of the 32 discovered); "Placement" = contact with SAM, substrate, both, or neither (distal),
                                with the minimum heavy-atom distances to the cofactor and to the substrate in angstroms
  motif_table.png               rendered table, banded by label
"""
import argparse
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CLASSES = ["DNA", "RNA", "Protein", "Metabolite"]
COLORS = {"DNA": "#6a51a3", "RNA": "#c14445", "Protein": "#618d43", "Metabolite": "#ed8c44"}

# Structural placement of each class-specific motif (structures/motif_placements.csv): minimum heavy-atom distance of the motif to the SAM/SAH
# cofactor and to the substrate in the best in-class experimental structure; "contact" means 4 A or less.
SPEC_MIN = 2.0       # class-specific: own-class carrier rate at least 2x the pooled rate outside the class (same cutoff as motif_structure_carriers.py)
PLACEMENT = {"near substrate/SAM": "SAM + Substrate", "near substrate": "Substrate", "near SAM": "SAM", "distal": "Distal"}


def parse_streme(cls, base=None):
    rows = []
    for m in re.finditer(r"MOTIF (\d+)-(\S+) STREME-\d+\s*\n"
                         r"letter-probability matrix: alength= \d+ w= (\d+) nsites= (\d+) E= (\S+)",
                         ((base or HERE) / f"streme_{cls}" / "streme.txt").read_text()):
        num, cons, w, nsites, ev = m.groups()
        rows.append({"num": int(num), "motif_id": f"{num}-{cons}", "consensus": cons,
                     "width": int(w), "streme_sites": int(nsites), "evalue": float(ev)})
    return rows


def fasta_count(path):
    return sum(1 for ln in open(path) if ln.startswith(">"))


def main():
    global HERE
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=None, help="<release>/motif_analysis (default: the archived table in this folder)")
    a = ap.parse_args()
    release = a.dir is not None
    base = Path(a.dir).resolve() if release else HERE
    HERE = base                                             # outputs go next to the inputs
    pl = pd.read_csv((base / "structure_search" / "motif_placements.csv") if release else (base / "structures" / "motif_placements.csv")).set_index("motif")
    place = pl["placement"].map(PLACEMENT)
    spec = pd.read_csv(base / "motif_specificity_matrix.csv").set_index("motif")
    fold = spec["ratio_vs_bg"]
    out = []
    for cls in CLASSES:
        if release:       # carriers of the motif's own class and their share of the class: the cross-scan's own-class column (same FIMO, same p < 1e-4)
            n_class = None
        else:
            n_class = fasta_count(base / f"class_{cls}.fasta")
            fimo = pd.read_csv(base / f"fimo_all_{cls}" / "fimo.tsv", sep="\t", comment="#")
            fimo = fimo.dropna(subset=["sequence_name"])
            carriers = fimo.groupby("motif_id")["sequence_name"].nunique()
        for r in parse_streme(cls, base):
            if release:
                n_prot = int(spec.loc[f"{cls}-{r['num']}", f"n_{cls}"]) if f"{cls}-{r['num']}" in spec.index else 0
                pct = float(spec.loc[f"{cls}-{r['num']}", f"pct_{cls}"]) if f"{cls}-{r['num']}" in spec.index else 0.0
            else:
                n_prot = int(carriers.get(r["motif_id"], 0))
                pct = round(100 * n_prot / n_class, 1)
            out.append({
                "Label": cls,
                "Motif": f"{cls}-{r['num']}",
                "Consensus": r["consensus"],
                "Width": r["width"],
                "E-value": f"{r['evalue']:.0e}",
                "Proteins": n_prot,
                "% of class": round(pct, 1),
                "Fold enrichment": f"{fold[f'{cls}-{r['num']}']:.1f}",
                "Placement": place.get(f"{cls}-{r['num']}", "—"),
                "To SAM (Å)": f"{pl.loc[f'{cls}-{r['num']}', 'd_SAM']:.1f}" if f"{cls}-{r['num']}" in pl.index else "—",
                "To substrate (Å)": f"{pl.loc[f'{cls}-{r['num']}', 'd_substrate']:.1f}" if f"{cls}-{r['num']}" in pl.index else "—",
                **({"Structure": pl.loc[f"{cls}-{r['num']}", "pdb"] if f"{cls}-{r['num']}" in pl.index else "—"} if release else {}),
            })
    df = pd.DataFrame(out)
    n_all = len(df)
    df = df[df["Fold enrichment"].astype(float) >= SPEC_MIN].reset_index(drop=True)      # drop the six pan-methyltransferase motifs below the cutoff
    print(f"[table] {len(df)} of {n_all} motifs reach {SPEC_MIN:g}-fold enrichment")
    df.to_csv(HERE / "motif_table.csv", index=False)
    print(df.to_string(index=False))
    print(f"\n[wrote] motif_table.csv ({len(df)} motifs)")

    # ---- rendered table ----
    cols = list(df.columns)
    fig, ax = plt.subplots(figsize=(15.0, 0.34 * len(df) + 0.3))
    ax.axis("off")
    tbl = ax.table(cellText=df.values, colLabels=cols, cellLoc="left",
                   loc="center", bbox=[0, 0, 1, 1])
    tbl.auto_set_font_size(False); tbl.set_fontsize(10)
    ncol = len(cols)
    widths = {"Label": 0.085, "Motif": 0.09, "Consensus": 0.135, "Width": 0.05, "E-value": 0.07, "Proteins": 0.07,
              "% of class": 0.08, "Fold enrichment": 0.11, "Placement": 0.11, "To SAM (Å)": 0.085, "To substrate (Å)": 0.115, "Structure": 0.08}
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_width(widths[cols[c]] if not release else widths[cols[c]] * 0.93)
        if r == 0:
            cell.set_facecolor("#333333"); cell.set_text_props(color="white", fontweight="bold")
        else:
            lab = df.iloc[r - 1]["Label"]
            if cols[c] == "Label":
                cell.set_facecolor(COLORS[lab]); cell.set_text_props(color="white", fontweight="bold")
            else:
                cell.set_facecolor("#ffffff" if (r % 2) else "#f4f4f4")
            if cols[c] == "Consensus":
                cell.set_text_props(fontfamily="monospace")
    fig.tight_layout()
    fig.savefig(HERE / "motif_table.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] motif_table.png")


if __name__ == "__main__":
    main()
