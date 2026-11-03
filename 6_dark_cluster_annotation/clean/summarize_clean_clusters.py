#!/usr/bin/env python3
"""Per dark protein-cluster CLEAN summary, from CLEAN's two raw per-protein outputs (run_clean_dark.sh).

  python summarize_clean_clusters.py --maxsep <name>_maxsep.csv --pvalue <name>_pvalue.csv --members <name>_members.csv --out <folder> [--ec-names ec_names.tsv | enzyme.dat]

  top_ec / top_ec_frac / n_distinct_ec   the most frequent first (best) maxsep call among the members, its share of the cluster, and the number of distinct first calls
  d_median                               the MEDIAN over the members of the distance of their best p-value call (CLEAN's distance to the EC cluster centre; lower = closer)
  confidence                             "reliable" d < 6, "marginal" 6 <= d < 8, "unreliable" d >= 8 (the cut points of the archived Tables S1/S2: with them all 33 archived
                                         clusters get the tier of the curated table, test_summarize_clean.py), written "tier (d=x.x)" as there
  mtase                                  "yes" if top_ec is an EC 2.1.1.* (a methyltransferase), else "No" (the archived curated column follows this rule for all 33 clusters)
  enzyme                                 the name of top_ec from --ec-names (a tab-separated EC, name table, e.g. from ExPASy ENZYME); empty without it
The curated table of the archived analysis also carries a biological-context note after an arrow ("-> anthracyclines"); that part is manual and is not produced here.
Writes clean_summary.csv to --out (columns: cluster, lineage, n, top_ec, top_ec_frac, n_distinct_ec, d_median, confidence, mtase, enzyme).
"""
import argparse
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

RELIABLE_BELOW, MARGINAL_BELOW = 6.0, 8.0


def read_calls(path):
    """CLEAN output rows 'Entry,EC:a.b.c.d/score[,EC:.../score]' -> {Entry: [(EC, score), ...]} in CLEAN's order (best first)."""
    out = {}
    for line in open(path):
        parts = line.strip().split(",")
        calls = [(m.group(1), float(m.group(2))) for x in parts[1:] if (m := re.match(r"EC:([\d.\-n]+)/([\d.eE+-]+)", x))]
        out[parts[0]] = calls
    return out


def tier(d):
    return "reliable" if d < RELIABLE_BELOW else "marginal" if d < MARGINAL_BELOW else "unreliable"


def summarize(maxsep, pvalue, members, ec_names=None):
    members = members.drop_duplicates("Entry")
    rows = []
    for c, g in members.groupby("cluster"):
        firsts = [maxsep[e][0][0] for e in g.Entry if maxsep.get(e)]
        ds = [pvalue[e][0][1] for e in g.Entry if pvalue.get(e)]
        if not firsts:
            continue
        top, n_top = Counter(firsts).most_common(1)[0]
        d = float(np.median(ds)) if ds else float("nan")
        rows.append({"cluster": int(c), "lineage": g.lineage.iloc[0] if "lineage" in g else "", "n": len(g), "top_ec": "EC:" + top, "top_ec_frac": round(n_top / len(firsts), 3),
                     "n_distinct_ec": len(set(firsts)), "d_median": round(d, 3), "confidence": f"{tier(d)} (d={d:.1f})", "mtase": "yes" if top.startswith("2.1.1.") else "No",
                     "enzyme": (ec_names or {}).get(top, "")})
    return pd.DataFrame(rows)


def read_enzyme_dat(path):
    """{EC: recommended name} from ExPASy's enzyme.dat (ID / DE lines per entry)."""
    names, ec, de = {}, None, []
    for line in open(path, errors="replace"):
        if line.startswith("ID   "):
            ec, de = line[5:].strip(), []
        elif line.startswith("DE   ") and ec:
            de.append(line[5:].strip())
        elif line.startswith("//") and ec:
            names[ec] = " ".join(de).rstrip(".")
            ec = None
    return names


def read_ec_names(path):
    """Either a tab-separated EC<TAB>name file or ExPASy's enzyme.dat."""
    if str(path).endswith(".dat"):
        return read_enzyme_dat(path)
    names = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            names[p[0].replace("EC:", "").strip()] = p[1].strip()
    return names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maxsep", required=True)
    ap.add_argument("--pvalue", required=True)
    ap.add_argument("--members", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ec-names")
    a = ap.parse_args()
    t = summarize(read_calls(a.maxsep), read_calls(a.pvalue), pd.read_csv(a.members), read_ec_names(a.ec_names) if a.ec_names else None)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "clean_summary.csv", index=False)
    print(f"{len(t)} clusters -> {out / 'clean_summary.csv'}")


if __name__ == "__main__":
    main()
