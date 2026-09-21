#!/usr/bin/env python3
"""Assemble the full STREME motif table: every discriminative motif for every
substrate label, with its consensus, discovery stats, and REAL per-protein
coverage (fraction of that class's members carrying it, from the FIMO-all run).

Inputs (this dir):
  streme_{cls}/streme.txt      8 STREME motifs per class (consensus, width, E-value)
  fimo_all_{cls}/fimo.tsv       FIMO hits of those motifs vs ALL class members
  class_{cls}.fasta             class member counts (denominator)
Outputs:
  motif_table.csv               one row per motif (32 rows)
  motif_table.png               rendered table, banded by label
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CLASSES = ["DNA", "RNA", "Protein", "Other"]
COLORS = {"DNA": "#6a51a3", "RNA": "#cb181d", "Protein": "#238b45", "Other": "#e6550d"}

# README-curated catalytic / cofactor roles, keyed by (class, motif number)
KNOWN = {
    ("DNA", 1): "C5-cytosine catalytic Pro-Cys (PCQ)",
    ("DNA", 7): "Motif IV [D/N]PPY (amino-MTase catalytic)",
    ("RNA", 1): "m5C catalytic Pro-Cys (DPPC)",
    ("RNA", 7): "basic RNA-binding patch (Lys-rich)",
    ("Protein", 5): "PRMT THW loop",
    ("Protein", 8): "NPPY (N-MTase catalytic)",
    ("Other", 3): "Motif I GxGxG (SAM-binding)",
}


def parse_streme(cls):
    rows = []
    for m in re.finditer(r"MOTIF (\d+)-(\S+) STREME-\d+\s*\n"
                         r"letter-probability matrix: alength= \d+ w= (\d+) nsites= (\d+) E= (\S+)",
                         (HERE / f"streme_{cls}" / "streme.txt").read_text()):
        num, cons, w, nsites, ev = m.groups()
        rows.append({"num": int(num), "motif_id": f"{num}-{cons}", "consensus": cons,
                     "width": int(w), "streme_sites": int(nsites), "evalue": float(ev)})
    return rows


def fasta_count(path):
    return sum(1 for ln in open(path) if ln.startswith(">"))


def main():
    out = []
    for cls in CLASSES:
        n_class = fasta_count(HERE / f"class_{cls}.fasta")
        fimo = pd.read_csv(HERE / f"fimo_all_{cls}" / "fimo.tsv", sep="\t", comment="#")
        fimo = fimo.dropna(subset=["sequence_name"])
        carriers = fimo.groupby("motif_id")["sequence_name"].nunique()
        for r in parse_streme(cls):
            n_prot = int(carriers.get(r["motif_id"], 0))
            out.append({
                "Label": cls,
                "Motif": f"{cls}-{r['num']}",
                "Consensus": r["consensus"],
                "Width": r["width"],
                "E-value": f"{r['evalue']:.0e}",
                "Proteins": n_prot,
                "% of class": round(100 * n_prot / n_class, 1),
                "Known role": KNOWN.get((cls, r["num"]), "—"),
            })
    df = pd.DataFrame(out)
    df.to_csv(HERE / "motif_table.csv", index=False)
    print(df.to_string(index=False))
    print(f"\n[wrote] motif_table.csv ({len(df)} motifs)")

    # ---- rendered table ----
    cols = list(df.columns)
    fig, ax = plt.subplots(figsize=(15, 0.34 * len(df) + 1.0))
    ax.axis("off")
    tbl = ax.table(cellText=df.values, colLabels=cols, cellLoc="left",
                   loc="center", bbox=[0, 0, 1, 1])
    tbl.auto_set_font_size(False); tbl.set_fontsize(10)
    ncol = len(cols)
    widths = {"Label": 0.07, "Motif": 0.075, "Consensus": 0.20, "Width": 0.05,
              "E-value": 0.08, "Proteins": 0.075, "% of class": 0.08, "Known role": 0.31}
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_width(widths[cols[c]])
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
    ax.set_title("Discriminative STREME motifs per substrate label — consensus, discovery E-value, "
                 "and real per-protein coverage (FIMO over all class members)",
                 fontsize=12, fontweight="bold", pad=14)
    fig.tight_layout()
    fig.savefig(HERE / "motif_table.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] motif_table.png")


if __name__ == "__main__":
    main()
