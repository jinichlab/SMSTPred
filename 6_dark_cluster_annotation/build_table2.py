#!/usr/bin/env python3
"""Build the manuscript's Tables S1 (Actinomycetota) and S2 (Plant) -- the two
independent computational calls for each of the 33 dark protein-clusters:
CLEAN's EC-centroid nearest-known-enzyme assignment, and Horizyn's dominant
reaction family from screening all 3,478 dark proteins against all 36,552
Rhea reactions (see ../README.md for the full methodology).

This reproduces the CSV content of Tables S1 and S2 as it appears in article_draft.docx
-- not the docx formatting/coloring itself, just the data.

Inputs:
  clean/clean_dark33_cluster_summary.csv  -- CLEAN side (manually curated: each
    cluster's dominant EC call -- see clean/clean_dark33_per_cluster_ec.csv for
    the raw automated call this was based on -- was looked up against enzyme
    nomenclature by hand to get a human-readable name and MTase yes/no judgment;
    that lookup step is not automated)
  horizyn/protein_clusters_family_final.csv  -- Horizyn side (fully automated)
  horizyn/avg_top1_score.csv  -- mean cosine-similarity score of each member's
    own rank-1 (top) reaction match, averaged per cluster -- a reliability
    signal independent of family-vote share (a cluster can vote unanimously
    for the same family while each member's best individual match is still a
    mediocre absolute score, or vice versa)
  family_descriptions.csv  -- short hand-written descriptions of each Horizyn
    reaction family's chemistry, verified against real Rhea equations

Output: table2_dark_clusters.csv

  python build_table2.py                       the archived inputs of this folder -> table2_dark_clusters.csv here (the 33 archived clusters)
  python build_table2.py --release-dir <dir>   a release's <release>/dark_annotation: clean/clean_summary.csv (clean/summarize_clean_clusters.py), horizyn/protein_clusters_family_final.csv
                                               + avg_top1_score.csv (horizyn/final_family_table.py) -> <dir>/table_dark_clusters.csv. The hand-written family descriptions
                                               (family_descriptions.csv) exist only for the families met so far; a family without one gets a DRAFT label made from the
                                               SAM-methylation flags of its Rhea reactions ("SAM methylation" / "NOT methylation") that has to be replaced by a described chemistry.
"""
import argparse
import csv
from pathlib import Path

import pandas as pd

IO = Path(__file__).parent.resolve()


def load_csv(path, key):
    with open(path) as f:
        return {row[key]: row for row in csv.DictReader(f)}


def draft_description(fam_id, rxn, titles=None):
    """No hand-written description: a label from the SAM-methylation flags of the family's Rhea reactions (rxn: rhea_reaction_clusters.csv), with one example reaction written
    from Rhea's equation text when `titles` ({rhea_lr_id: 'substrate -> product'}) is given."""
    g = rxn[rxn.reaction_cluster == int(fam_id)]
    k = int(g.is_sam_methylation.sum())
    kind = "SAM methylation" if k >= len(g) / 2 else "NOT methylation"
    ex = next((titles[i] for i in g.sort_values(["is_sam_methylation", "rhea_lr_id"], ascending=[False, True]).rhea_lr_id if titles and titles.get(i)), None)
    return f"[draft] {kind}" + (f", e.g. {ex}" if ex else "") + f" ({k} of {len(g)} reactions are SAM methylation)"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release-dir")
    ap.add_argument("--descriptions", default=str(IO / "family_descriptions.csv"))
    ap.add_argument("--rxn-clusters", default=str(IO / "horizyn" / "rhea_reaction_clusters.csv"))
    a = ap.parse_args()
    if a.release_dir:
        rd = Path(a.release_dir)
        clean_path, fam_path, avg_path, out = rd / "clean" / "clean_summary.csv", rd / "horizyn" / "protein_clusters_family_final.csv", rd / "horizyn" / "avg_top1_score.csv", rd / "table_dark_clusters.csv"
    else:
        clean_path, fam_path, avg_path, out = (IO / "clean" / "clean_dark33_cluster_summary.csv", IO / "horizyn" / "protein_clusters_family_final.csv",
                                               IO / "horizyn" / "avg_top1_score.csv", IO / "table2_dark_clusters.csv")
    clean = load_csv(clean_path, "cluster")
    fam = load_csv(fam_path, "cluster")
    desc = load_csv(a.descriptions, "family_id")
    avg_score = load_csv(avg_path, "cluster")
    rxn = pd.read_csv(a.rxn_clusters) if a.release_dir else None
    titles, notes = None, {}
    if a.release_dir:
        rel = rd.parent                                    # <release>/dark_annotation -> <release>
        eqf = rel / "reference_tables" / "rhea_equations.tsv"
        if eqf.exists():                                   # one example reaction per family, written from Rhea's equation text
            from rhea_text import load_equations, reaction_title
            eq = load_equations(eqf)
            dd = pd.read_csv(rel / "reference_tables" / "rhea-directions.tsv", sep="\t")
            master = {}
            for col in ("RHEA_ID_MASTER", "RHEA_ID_LR", "RHEA_ID_RL", "RHEA_ID_BI"):
                master.update(zip(dd[col].astype(int), dd["RHEA_ID_MASTER"].astype(int)))
            titles = {int(i): reaction_title(eq.get(master.get(int(i), int(i)), "")) for i in rxn.rhea_lr_id}
        cn = rd / "context_notes.csv"
        if cn.exists():                                    # the automatic pathway note of the nearest characterized enzyme (context_notes.py)
            notes = {int(r.cluster): r.pathway for r in pd.read_csv(cn).fillna("").itertuples() if r.pathway}

    rows = []
    for cluster, c in clean.items():
        f = fam.get(cluster)
        if f is None:
            if cluster in avg_score:   # scored, but every member's rank-1 reaction is HDBSCAN noise: no family to report
                horizyn_family, horizyn_share = "(no family: all top-1 reactions are noise)", "0% of cluster, 100% noise"
            else:
                horizyn_family, horizyn_share = "(not in Horizyn screen)", "-"
        else:
            fam_id = f["coherent_dominant_family"]
            d = desc.get(fam_id) or {"description": draft_description(fam_id, rxn, titles) if rxn is not None else "(see repository)"}
            horizyn_family = (f"Family {fam_id}: {d['description']} "
                               f"({f['coherent_family_size']} rxns, "
                               f"{float(f['coherent_family_pct_orphan']):.0f}% orphan)")
            horizyn_share = (f"{float(f['coherent_family_share']) * 100:.0f}% of cluster, "
                              f"{float(f['pct_noise_votes']):.0f}% noise")
        a = avg_score.get(cluster)
        horizyn_avg_score = f"{float(a['avg_top1_score']):.2f}" if a else "-"
        rows.append(dict(
            cluster=cluster, lineage=c["lineage"], n=c["n"],
            clean_enzyme=(c["enzyme"] + f"  → {notes[int(cluster)]}" if int(cluster) in notes and not str(c["confidence"]).startswith("unreliable") else c["enzyme"]), clean_mtase=c["mtase"], clean_confidence=c["confidence"],
            horizyn_family=horizyn_family, horizyn_share=horizyn_share,
            horizyn_avg_score=horizyn_avg_score,
        ))

    rows.sort(key=lambda r: (r["lineage"], -int(r["n"])))
    fields = ["cluster", "lineage", "n", "clean_enzyme", "clean_mtase", "clean_confidence",
              "horizyn_family", "horizyn_share", "horizyn_avg_score"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} clusters -> {out}")


if __name__ == "__main__":
    main()
