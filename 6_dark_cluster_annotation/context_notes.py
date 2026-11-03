#!/usr/bin/env python3
"""Automatic 'biological context' note for every dark cluster (what used to be typed into Tables S1 / S2, e.g. "-> anthracyclines").

  python context_notes.py --release <release>

For the EC number CLEAN assigned to the cluster (clean/clean_summary.csv: top_ec) the reviewed (Swiss-Prot) UniProt entries with that EC are asked for their PATHWAY annotation
(rest.uniprot.org; answers cached in dark_annotation/context_cache.json); the most frequent pathway, with the number of reviewed enzymes that carry it, is the note. If the EC has no
pathway annotation the note is empty. Horizyn's most voted reaction of the cluster, written as 'substrate -> product' from Rhea's equation text, is added as a second, independent hint.
Writes <release>/dark_annotation/context_notes.csv: cluster, ec, n_reviewed_with_ec, pathway, pathway_share, horizyn_top_reaction, note
build_table2.py appends the pathway to the enzyme name of the table ("aklanonic acid methyltransferase -> Antibiotic biosynthesis; anthracycline biosynthesis").
An automatic note is a lead, not a finding: it says where the closest characterized enzyme works, which holds only as far as CLEAN's assignment does (see the confidence column).
"""
import argparse
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from rhea_text import load_equations, reaction_title  # noqa: E402


def pathways_of(ec, cache):
    """[(pathway text, n reviewed entries)] for a reviewed-entry search by EC; cached."""
    if ec in cache:
        return cache[ec]
    q = f"(ec:{ec}) AND (reviewed:true)"
    for attempt in range(4):
        try:
            r = requests.get("https://rest.uniprot.org/uniprotkb/search", params={"query": q, "fields": "accession,cc_pathway", "format": "tsv", "size": 500}, timeout=120)
            r.raise_for_status()
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))
    lines = r.text.splitlines()[1:]
    c = Counter()
    for ln in lines:
        txt = ln.split("\t", 1)[1] if "\t" in ln else ""
        for p in txt.split("PATHWAY:"):                    # one entry can list several PATHWAY comments
            p = re.sub(r"\s*\{ECO:[^}]*\}", "", p)
            p = re.sub(r"\.\.+", ".", p)                     # UniProt writes "..;" where a comment ends in a full stop
            p = re.sub(r"[.;\s]+$", "", p.strip())
            if p:
                c[p] += 1
    cache[ec] = {"n_reviewed": len(lines), "pathways": c.most_common(5)}
    return cache[ec]


def horizyn_top(topk, members, eq, directions):
    """{cluster: 'substrate -> product' of the reaction most members voted for (rank 1)}"""
    t = topk[topk["rank"] == 1].merge(members[["Entry", "cluster"]], left_on="accession", right_on="Entry", suffixes=("_x", ""))
    out = {}
    for c, g in t.groupby("cluster"):
        rid = g.rhea_lr_id.value_counts().index[0]
        out[c] = reaction_title(eq.get(directions.get(int(rid), int(rid)), "")) or ""
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", required=True)
    a = ap.parse_args()
    rel = Path(a.release)
    A = rel / "dark_annotation"
    name = (A / "dark_name.txt").read_text().strip()
    clean = pd.read_csv(A / "clean" / "clean_summary.csv")
    members = pd.read_csv(A / f"{name}_members.csv")
    topk = pd.read_csv(A / "horizyn_run" / "horizyn_full_topk.csv")
    eq = load_equations(rel / "reference_tables" / "rhea_equations.tsv")
    d = pd.read_csv(rel / "reference_tables" / "rhea-directions.tsv", sep="\t")
    master = {}
    for col in ("RHEA_ID_MASTER", "RHEA_ID_LR", "RHEA_ID_RL", "RHEA_ID_BI"):
        master.update(zip(d[col].astype(int), d["RHEA_ID_MASTER"].astype(int)))
    top = horizyn_top(topk, members, eq, master)
    cf = A / "context_cache.json"
    cache = json.load(open(cf)) if cf.exists() else {}
    rows = []
    for r in clean.itertuples():
        ec = str(r.top_ec).replace("EC:", "")
        info = pathways_of(ec, cache)
        pw, n = (info["pathways"][0] if info["pathways"] else ("", 0))
        rows.append({"cluster": r.cluster, "ec": ec, "n_reviewed_with_ec": info["n_reviewed"], "pathway": pw, "pathway_share": f"{n}/{info['n_reviewed']}" if info["n_reviewed"] else "",
                     "horizyn_top_reaction": top.get(r.cluster, ""), "note": pw})
    json.dump(cache, open(cf, "w"))
    out = pd.DataFrame(rows)
    out.to_csv(A / "context_notes.csv", index=False)
    print(out[["cluster", "ec", "n_reviewed_with_ec", "pathway", "pathway_share", "horizyn_top_reaction"]].to_string(index=False, max_colwidth=60))
    print(f"\n{(out.pathway != '').sum()} of {len(out)} clusters have a pathway note -> {A / 'context_notes.csv'}")


if __name__ == "__main__":
    main()
