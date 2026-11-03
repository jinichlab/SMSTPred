"""Which dark clusters Fig 8 shows, chosen by a rule instead of by hand.

A cluster is shown when CLEAN and Horizyn agree on its chemistry at the level of the reaction family:
  * CLEAN is confident: tier "reliable" (median best-call distance d < 6) and the call is a methyltransferase (EC 2.1.1.*), and
  * Horizyn's dominant reaction family for the cluster (final_family_table.py) contains a reaction that carries CLEAN's EC (Rhea's rhea2ec.tsv).
Clusters are listed by number; at most `max_panels` (six) are drawn, the largest first when there are more.
On release 2026_03 this picks c317, c353, c512 and c729. The rule is deliberately generous at the family level (a family can mix several ECs):
the figure then shows what Horizyn actually voted for, which is the evidence the reader weighs.
"""
import pandas as pd


def master_map(directions):
    """{any Rhea entry id: its master id} from rhea-directions.tsv (master, left->right, right->left, bidirectional). The master is NOT a fixed offset from its entries (15573 -> 15574, 15575, 15576)."""
    m = {}
    for col in ("RHEA_ID_MASTER", "RHEA_ID_LR", "RHEA_ID_RL", "RHEA_ID_BI"):
        m.update(zip(directions[col].astype(int), directions["RHEA_ID_MASTER"].astype(int)))
    return m


def family_ecs(rxn_clusters, rhea2ec, masters):
    """{reaction family: set of ECs of its reactions}."""
    ec = rhea2ec.groupby("MASTER_ID")["ID"].apply(set).to_dict()
    r = rxn_clusters[rxn_clusters["reaction_cluster"] != -1]
    out = {}
    for fam, rid in zip(r["reaction_cluster"], r["rhea_lr_id"]):
        out.setdefault(int(fam), set()).update(ec.get(str(masters.get(int(rid), "")), set()))
    return out


def select(clean, fam, rxn_clusters, rhea2ec, masters, max_panels=6):
    """clean: clean_summary.csv; fam: protein_clusters_family_final.csv. Returns the selected cluster numbers (list of int)."""
    fe = family_ecs(rxn_clusters, rhea2ec, masters)
    d = clean.merge(fam[["cluster", "coherent_dominant_family"]], on="cluster", how="inner")
    keep = []
    for _, r in d.iterrows():
        if not str(r["confidence"]).startswith("reliable") or r["mtase"] != "yes":
            continue
        if str(r["top_ec"]).replace("EC:", "") in fe.get(int(r["coherent_dominant_family"]), set()):
            keep.append((int(r["cluster"]), int(r["n"])))
    if len(keep) > max_panels:
        keep = sorted(keep, key=lambda t: -t[1])[:max_panels]
    return sorted(c for c, _ in keep)
