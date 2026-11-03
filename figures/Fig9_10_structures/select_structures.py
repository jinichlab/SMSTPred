"""Which motif structures Figs 9 and 10 show, chosen by rule from the release's motif placements (7_motif_analysis/motif_structure_pick.py).

Fig 9 (known sites): four motifs whose structure has the motif in contact with both SAM and the substrate ("near substrate/SAM"), the CLOSEST of each class
  (closeness = the larger of the two distances, motif -> SAM and motif -> substrate; ties: their sum).
Fig 10 (distal sites): four motifs in contact with neither ("distal"), the FURTHEST of each class (remoteness = the smaller of the two distances; ties: their sum).
One per class where the class has one (order: DNA, RNA, Protein, Metabolite); a slot that a class cannot fill is filled with the next closest (Fig 9) / furthest (Fig 10) motif of ANY class
that is not drawn yet. If fewer than four motifs have the placement at all, the remaining slots are filled with the closest / furthest of the others (flagged in `note`).
Motifs whose structure is a weak pick (a carrier of another class, or a substrate that is not of the class's own kind: `in_class` False or `substrate_kind_fallback` True) are not used unless
allow_weak is set: a figure should not rest on them.
"""
import pandas as pd

CLASSES = ["DNA", "RNA", "Protein", "Metabolite"]


def closeness(r):
    return (max(r["d_SAM"], r["d_substrate"]), r["d_SAM"] + r["d_substrate"])


def remoteness(r):
    return (-min(r["d_SAM"], r["d_substrate"]), -(r["d_SAM"] + r["d_substrate"]))


def choose(pl, wanted_placement, key, n=4, allow_weak=False):
    """pl: motif_placements.csv. Returns a list of dicts (motif, why) of up to n motifs."""
    pl = pl.copy()
    weak = (~pl["in_class"].astype(bool)) | pl["substrate_kind_fallback"].astype(bool)
    if not allow_weak:
        pl = pl[~weak]
    picked, why = [], {}
    for cls in CLASSES:
        c = pl[(pl["class"] == cls) & (pl["placement"] == wanted_placement)]
        if len(c):
            m = min(c.to_dict("records"), key=key)
            picked.append(m["motif"]); why[m["motif"]] = f"best {cls} motif"
    pool = pl[~pl["motif"].isin(picked)]
    for stage, sub in (("next overall", pool[pool["placement"] == wanted_placement]), ("next overall, other placement", pool[pool["placement"] != wanted_placement])):
        while len(picked) < n and len(sub):
            m = min(sub.to_dict("records"), key=key)
            picked.append(m["motif"]); why[m["motif"]] = stage
            sub = sub[sub["motif"] != m["motif"]]
    order = {c: i for i, c in enumerate(CLASSES)}
    out = sorted(picked, key=lambda m: (order[m.split("-")[0]], int(m.split("-")[1])))
    return [{"motif": m, "why": why[m]} for m in out]


def select(pl, allow_weak=False):
    known = choose(pl, "near substrate/SAM", closeness, 4, allow_weak)
    distal = choose(pl, "distal", remoteness, 4, allow_weak)
    for r in known:
        r["figure"], r["panel_set"] = 9, "Known"
    for r in distal:
        r["figure"], r["panel_set"] = 10, "Distal"
    sel = pd.DataFrame(known + distal)
    return sel.merge(pl, on="motif", how="left")


def caption(r):
    """Enzyme name (UniProt's, without the EC number and the synonym brackets), the motif and its matched sequence, the PDB entry."""
    import re
    name = re.sub(r"\s*\(EC [^)]*\).*$", "", str(r["protein"]))          # UniProt lists "Name (EC x) (synonym) ...": keep the name
    name = name if len(name) <= 46 else name[:44] + "…"
    return f"{name} · {r['motif']} {r['matched_seq']} ({r['pdb']})"
