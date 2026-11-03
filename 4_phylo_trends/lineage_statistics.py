#!/usr/bin/env python3
"""Lineage comparisons quoted in the Results (plants and Actinomycetota as small-molecule specialists; Micrococcales; Chlorophyta).

Per-organism SAM-methyltransferase counts by predicted substrate come from figures/Fig5_proteome_scaling/proteome_vs_sam_bylabel.csv; the
lineage of each organism from data/selected_1259_organisms.csv. Two-sided Mann-Whitney U tests, not corrected for multiple comparisons.

Outputs (next to this script): lineage_statistics.csv (one row per comparison and measure) and lineage_statistics.txt (the numbers as quoted).

  python lineage_statistics.py                    the archived table (figures/Fig5_proteome_scaling/proteome_vs_sam_bylabel.csv), output next to this script
  python lineage_statistics.py --release <dir>    a refresh release: <dir>/fig5/proteome_vs_sam_bylabel.csv (figures/Fig5_proteome_scaling/run_fig5.sh), output into <dir>/fig5/
"""
import argparse
import os
import re
from pathlib import Path

import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = Path(os.environ.get("SMSTPRED_DATA", REPO / "data"))


def canon(name):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9_\-]+", "", re.sub(r"\s+", "_", str(name).strip()))).strip("_")


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--release", default=None, help="release folder (default: the archived table)")
args = ap.parse_args()
if args.release:
    SRC, OUT = Path(args.release) / "fig5" / "proteome_vs_sam_bylabel.csv", Path(args.release) / "fig5"
else:
    SRC, OUT = REPO / "figures/Fig5_proteome_scaling/proteome_vs_sam_bylabel.csv", HERE
per = pd.read_csv(SRC)
meta = pd.read_csv(DATA / "selected_1259_organisms.csv")
per["key"], meta["key"] = per.organism.map(canon), meta.Organism.map(canon)
d = per.merge(meta[["key", "Taxonomic lineage"]].drop_duplicates("key"), on="key", how="left")
print(f"{d['Taxonomic lineage'].notna().sum()} of {len(d)} organisms have a lineage in selected_1259_organisms.csv", flush=True)
lin = d["Taxonomic lineage"].fillna("")
d["pct_of_proteome"] = 100 * d.n_sam / d.proteome_size
d["frac_metabolite"] = 100 * d.n_Metabolite / d.n_sam
d = d[d.n_sam > 0]
lin = d["Taxonomic lineage"].fillna("")
G = {
    "Micrococcales": d[lin.str.contains("Micrococcales")],
    "other Actinomycetota": d[lin.str.contains("Actinomycetota") & ~lin.str.contains("Micrococcales")],
    "Actinomycetota": d[lin.str.contains("Actinomycetota")],
    "other bacteria": d[(d.domain == "Bacteria") & ~lin.str.contains("Actinomycetota")],
    "Chlorophyta": d[lin.str.contains("Chlorophyta")],
    "vascular plants": d[lin.str.contains("Tracheophyta")],
    "plants": d[lin.str.contains("Viridiplantae")],
    "other eukaryotes": d[(d.domain == "Eukaryota") & ~lin.str.contains("Viridiplantae")],
}
COMPARE = [("Micrococcales", "other Actinomycetota"), ("Chlorophyta", "vascular plants"), ("plants", "other eukaryotes"),
           ("Actinomycetota", "other bacteria")]
MEASURES = ["n_sam", "n_Metabolite", "n_RNA", "n_Protein", "n_DNA", "frac_metabolite", "pct_of_proteome"]
rows = []
for a, b in COMPARE:
    for m in MEASURES:
        x, y = G[a][m], G[b][m]
        rows.append(dict(group_a=a, n_a=len(x), median_a=x.median(), group_b=b, n_b=len(y), median_b=y.median(),
                         measure=m, mannwhitney_U=stats.mannwhitneyu(x, y).statistic, p=stats.mannwhitneyu(x, y).pvalue))
res = pd.DataFrame(rows)
res.to_csv(OUT / "lineage_statistics.csv", index=False, float_format="%.6g")
q = res.set_index(["group_a", "measure"])
def line(a, m, label, fmt="{:.1f}"):
    r = q.loc[(a, m)]
    return f"{label}: {a} (n={r.n_a}) median {fmt.format(r.median_a)} vs {r.group_b} (n={r.n_b}) median {fmt.format(r.median_b)}, p = {r.p:.1e}"
lines = [line("Micrococcales", "n_sam", "SAM-MTases per organism"), line("Micrococcales", "frac_metabolite", "Metabolite fraction (%)"),
         line("Micrococcales", "n_Metabolite", "Metabolite count"), line("Micrococcales", "n_RNA", "RNA count"),
         line("Chlorophyta", "n_sam", "SAM-MTases per organism"), line("Chlorophyta", "frac_metabolite", "Metabolite fraction (%)"),
         line("Chlorophyta", "pct_of_proteome", "SAM-MTases, % of proteome", "{:.2f}"),
         line("plants", "n_sam", "SAM-MTases per organism"), line("plants", "n_Metabolite", "Metabolite count"), line("plants", "frac_metabolite", "Metabolite fraction (%)"),
         line("Actinomycetota", "frac_metabolite", "Metabolite fraction (%)")]
(OUT / "lineage_statistics.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
