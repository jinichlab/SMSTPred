#!/usr/bin/env python3
"""Aggregate CLEAN's raw per-protein EC calls (dark33_maxsep.csv, from
run_clean_dark33.sh) up to the 33 dark protein-clusters, for transparency /
cross-checking against the manually-curated `clean_dark33_cluster_summary.csv`
(which additionally assigns a human-readable enzyme name and an MTase yes/no
judgment per cluster -- that step is manual EC-nomenclature lookup, not
automated here).

Input:  dark33_maxsep.csv (Entry, EC:call[/score][,EC:call[/score]...] --
          CLEAN's maxsep-thresholded top call(s) per protein), ../dark33_members.csv
Output: clean_dark33_per_cluster_ec.csv (cluster, n, top_ec, top_ec_frac,
          n_distinct_ec)
"""
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

IO = Path(__file__).parent.resolve()


def load_top_ec(path):
    """First EC:call in each row is CLEAN's #1 (maxsep) pick for that protein."""
    top_ec = {}
    with open(path) as f:
        for line in f:
            entry, rest = line.strip().split(",", 1)
            first_call = rest.split(",")[0]
            m = re.match(r"(EC:[\d.]+)", first_call)
            top_ec[entry] = m.group(1) if m else first_call
    return top_ec


def main():
    top_ec = load_top_ec(IO / "dark33_maxsep.csv")

    cluster_of = {}
    with open(IO.parent / "dark33_members.csv") as f:
        for row in csv.DictReader(f):
            cluster_of[row["Entry"]] = int(row["cluster"])

    by_cluster = defaultdict(list)
    for entry, ec in top_ec.items():
        c = cluster_of.get(entry)
        if c is not None:
            by_cluster[c].append(ec)

    rows = []
    for c, ecs in sorted(by_cluster.items()):
        counts = Counter(ecs)
        top, n_top = counts.most_common(1)[0]
        rows.append(dict(cluster=c, n=len(ecs), top_ec=top,
                          top_ec_frac=round(n_top / len(ecs), 3),
                          n_distinct_ec=len(counts)))

    out = IO / "clean_dark33_per_cluster_ec.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cluster", "n", "top_ec", "top_ec_frac", "n_distinct_ec"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} clusters -> {out.name}")


if __name__ == "__main__":
    main()
