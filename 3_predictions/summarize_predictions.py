#!/usr/bin/env python3
"""Summary statistics of the tree-of-life predictions quoted in the manuscript (Results, Stage 3).

Reads the per-organism CSVs in predictions_output/ and the training labels, and writes prediction_summary.txt / .csv:
how many proteins were scored, how many were also in the training set, how many received exactly one class, several classes
or none, how many fall outside the applicability domain, and the ensemble-disagreement medians.
"""
from pathlib import Path
import os

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("SMSTPRED_DATA", HERE.parent / "data"))
CLASSES = ["Metabolite", "RNA", "Protein", "DNA"]

frames = []
for f in sorted((HERE / "predictions_output").glob("*_predictions.csv")):
    d = pd.read_csv(f)
    d["organism_file"] = f.name
    frames.append(d)
p = pd.concat(frames, ignore_index=True)
train = set(pd.read_csv(DATA / "labels_SAM_full.csv", usecols=["Entry"]).Entry)

calls = p[[f"call_{c}" for c in CLASSES]].to_numpy().astype(int)
n_calls = calls.sum(axis=1)
p["n_calls"] = n_calls
p["in_training"] = p.Entry.isin(train)
n = len(p)

out = []
def say(label, value, note=""):
    out.append((label, value, note))

say("organism files", len(frames))
say("proteins scored (rows)", n)
say("unique accessions", p.Entry.nunique())
say("rows whose accession is in the SAM training set", int(p.in_training.sum()), "labelled training proteins that are also in the proteomes")
say("rows not in the training set", int((~p.in_training).sum()))
for c in CLASSES:
    k = int(((n_calls == 1) & (p[f"call_{c}"] == 1)).sum())
    say(f"exactly one call: {c}", k, f"{100 * k / n:.1f}% of rows")
k = int((n_calls > 1).sum())
say("two or more calls (multi-label)", k, f"{100 * k / n:.1f}% of rows")
k = int((n_calls == 0).sum())
say("no call (every score below its threshold)", k, f"{100 * k / n:.1f}% of rows")
k = int((p.in_domain == 0).sum())
say("outside the applicability domain (in_domain = 0)", k, f"{100 * k / n:.1f}% of rows")
for name, mask in [("all", np.ones(n, bool)), ("in-domain", (p.in_domain == 1).to_numpy()), ("out-of-domain", (p.in_domain == 0).to_numpy())]:
    say(f"median ensemble disagreement, {name}", round(float(p.disagreement[mask].median()), 4))
x = p[(p.in_training == False)]
say("median ad_maxsim, proteins not in the training set", round(float(x.ad_maxsim.median()), 3))

res = pd.DataFrame(out, columns=["quantity", "value", "note"])
res.to_csv(HERE / "prediction_summary.csv", index=False)
width = max(len(a) for a, _, _ in out)
lines = [f"{a:<{width}s}  {b:>10}  {c}" for a, b, c in out]
(HERE / "prediction_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
