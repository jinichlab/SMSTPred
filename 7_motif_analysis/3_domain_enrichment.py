#!/usr/bin/env python3
"""
Stage 6 motif analysis -- step 3 (A): InterPro domain enrichment per class.

For each single-label substrate class, find InterPro domains over-represented
vs all other classes (one-sided Fisher exact on presence/absence, plus
fold-enrichment). Grounds "what distinguishes the classes" at the domain level.

  python 3_domain_enrichment.py                 reads/writes this folder (the archived corpus.csv)
  python 3_domain_enrichment.py --dir <folder>  reads <folder>/corpus.csv, writes <folder>/domain_enrichment.csv (a release: <release>/motif_analysis,
                                                made by 1_build_corpus.py --release)

Output: domain_enrichment.csv  (class, interpro, name, n_in, pct_in, pct_out,
                                 fold, p), and a printed top-12 per class.
"""
import argparse
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

IO = Path(__file__).resolve().parent
LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
MIN_IN = 20   # minimum occurrences in-class to test a domain


def parse_interpro(field):
    """Return {IPR id: name} from an InterPro cell like
    'IPR029063 SAM-dependent MTase;IPR...' (names optional)."""
    out = {}
    for tok in str(field).split(";"):
        m = re.match(r"\s*(IPR\d+)\s*(.*)", tok)
        if m:
            out[m.group(1)] = m.group(2).strip()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=str(IO), help="folder with corpus.csv; domain_enrichment.csv is written there")
    ap.add_argument("--entry-list", default=str(IO / "interpro_entry.list"), help="InterPro entry.list (accession, type, name)")
    a = ap.parse_args()
    out = Path(a.dir)
    d = pd.read_csv(out / "corpus.csv", low_memory=False)
    d = d[d["class"].isin(LABELS)].copy()
    d["ipr"] = d["InterPro"].fillna("").apply(parse_interpro)
    print("sample InterPro:", d["InterPro"].dropna().iloc[0][:90], flush=True)

    name_map = {}
    with open(a.entry_list) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                name_map[p[0]] = p[2]
    N = len(d)
    all_dom = set().union(*d["ipr"].apply(lambda m: set(m))) if N else set()

    rows = []
    for lab in LABELS:
        inmask = (d["class"] == lab).values
        nin, nout = int(inmask.sum()), int((~inmask).sum())
        cin, cout = Counter(), Counter()
        for ipr, isin in zip(d["ipr"], inmask):
            (cin if isin else cout).update(ipr.keys())
        for dom in all_dom:
            a, b = cin[dom], cout[dom]
            if a < MIN_IN:
                continue
            fin, fout = a / nin, (b / nout if nout else 0.0)
            fold = fin / fout if fout > 0 else np.inf
            _, p = fisher_exact([[a, nin - a], [b, nout - b]], alternative="greater")
            rows.append({"class": lab, "interpro": dom, "name": name_map.get(dom, ""),
                         "n_in": a, "pct_in": round(fin, 4), "pct_out": round(fout, 4),
                         "fold": round(fold, 2), "p": p})
    m = pd.DataFrame(rows).sort_values(["class", "fold", "n_in", "interpro"], ascending=[True, False, False, False], kind="stable")  # ties: InterPro accession, descending
    m.to_csv(out / "domain_enrichment.csv", index=False)

    for lab in LABELS:
        print(f"\n=== {lab}: most prevalent class-specific InterPro domains ===", flush=True)
        sub = m[(m["class"] == lab) & (m["fold"] >= 2)].sort_values("n_in", ascending=False).head(12)
        for _, r in sub.iterrows():
            fold = "inf" if np.isinf(r.fold) else f"{r.fold:.1f}x"
            print(f"  {r.interpro}  {fold:>6}  in={r.pct_in:.0%} out={r.pct_out:.0%}  "
                  f"n={int(r.n_in)}  {r['name'][:48]}", flush=True)
    print("\nwrote domain_enrichment.csv", flush=True)


if __name__ == "__main__":
    main()
