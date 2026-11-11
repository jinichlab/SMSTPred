#!/usr/bin/env python3
"""Out-of-fold accuracy of the deployed classifier on enzymes with and without a class-specific STREME motif.

  python motif_free_accuracy.py --release <release>      (reads deployment/oof_predictions.csv and motif_analysis/motif_status_all.csv)

Single-label benchmark proteins only; 'has_motif' = the protein carries at least one FIMO hit of a motif of its own class (fimo_all_manifest.py).
Writes <release>/motif_analysis/motif_free_accuracy.csv (per class and overall: n, accuracy for has_motif True/False).
"""
import argparse
from pathlib import Path
import pandas as pd

L = ["Metabolite", "Protein", "RNA", "DNA"]
ap = argparse.ArgumentParser(); ap.add_argument("--release", required=True); a = ap.parse_args()
r = Path(a.release)
o = pd.read_csv(r / "deployment/oof_predictions.csv")
s = pd.read_csv(r / "motif_analysis/motif_status_all.csv")
o = o[o[L].sum(axis=1) == 1].copy()
o["true"] = o[L].idxmax(axis=1)
o["pred"] = [L[i] for i in o[["oof_" + k for k in L]].values.argmax(1)]
o["ok"] = o.true == o.pred
m = o.merge(s[["Entry", "class", "has_motif"]], on="Entry")
m = m[m["class"] == m["true"]]
rows = []
for k, g in [("all", m)] + [(c, m[m.true == c]) for c in L]:
    for h, gg in g.groupby("has_motif"):
        rows.append(dict(cls=k, has_motif=h, n=len(gg), accuracy=round(gg.ok.mean(), 4)))
out = pd.DataFrame(rows); out.to_csv(r / "motif_analysis/motif_free_accuracy.csv", index=False); print(out.to_string(index=False))
