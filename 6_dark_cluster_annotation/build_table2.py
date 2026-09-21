#!/usr/bin/env python3
"""Build the manuscript's Table 2 (a: Actinomycetota, b: Plant) -- the two
independent computational calls for each of the 33 dark protein-clusters:
CLEAN's EC-centroid nearest-known-enzyme assignment, and Horizyn's dominant
reaction family from screening all 3,478 dark proteins against all 36,552
Rhea reactions (see ../README.md for the full methodology).

This reproduces the CSV content of Table 2 as it appears in article_draft.docx
-- not the docx formatting/coloring itself, just the data.

Inputs:
  clean/clean_dark33_cluster_summary.csv  -- CLEAN side (manually curated: each
    cluster's dominant EC call -- see clean/clean_dark33_per_cluster_ec.csv for
    the raw automated call this was based on -- was looked up against enzyme
    nomenclature by hand to get a human-readable name and MTase yes/no judgment;
    that lookup step is not automated)
  horizyn/protein_clusters_family_final.csv  -- Horizyn side (fully automated)
  family_descriptions.csv  -- short hand-written descriptions of each Horizyn
    reaction family's chemistry, verified against real Rhea equations

Output: table2_dark_clusters.csv
"""
import csv
from pathlib import Path

IO = Path(__file__).parent.resolve()


def load_csv(path, key):
    with open(path) as f:
        return {row[key]: row for row in csv.DictReader(f)}


def main():
    clean = load_csv(IO / "clean" / "clean_dark33_cluster_summary.csv", "cluster")
    fam = load_csv(IO / "horizyn" / "protein_clusters_family_final.csv", "cluster")
    desc = load_csv(IO / "family_descriptions.csv", "family_id")

    rows = []
    for cluster, c in clean.items():
        f = fam.get(cluster)
        if f is None:
            horizyn_family, horizyn_share = "(not in Horizyn screen)", "-"
        else:
            fam_id = f["coherent_dominant_family"]
            d = desc.get(fam_id, {"description": "(see repository)"})
            horizyn_family = (f"Family {fam_id}: {d['description']} "
                               f"({f['coherent_family_size']} rxns, "
                               f"{float(f['coherent_family_pct_orphan']):.0f}% orphan)")
            horizyn_share = (f"{float(f['coherent_family_share']) * 100:.0f}% of cluster, "
                              f"{float(f['pct_noise_votes']):.0f}% noise")
        rows.append(dict(
            cluster=cluster, lineage=c["lineage"], n=c["n"],
            clean_enzyme=c["enzyme"], clean_mtase=c["mtase"], clean_confidence=c["confidence"],
            horizyn_family=horizyn_family, horizyn_share=horizyn_share,
        ))

    rows.sort(key=lambda r: (r["lineage"], -int(r["n"])))
    fields = ["cluster", "lineage", "n", "clean_enzyme", "clean_mtase", "clean_confidence",
              "horizyn_family", "horizyn_share"]
    out = IO / "table2_dark_clusters.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} clusters -> {out.name}")


if __name__ == "__main__":
    main()
