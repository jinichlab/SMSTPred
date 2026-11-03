#!/usr/bin/env python3
"""select_structures.choose / select on a made-up placement table."""
import pandas as pd
from select_structures import choose, closeness, remoteness, select


def row(motif, place, dc, ds, in_class=True, fb=False):
    return {"motif": motif, "class": motif.split("-")[0], "placement": place, "d_SAM": dc, "d_substrate": ds, "in_class": in_class, "substrate_kind_fallback": fb,
            "protein": "P", "pdb": "1ABC", "matched_seq": "AAA"}


S = "near substrate/SAM"
pl = pd.DataFrame([row("DNA-1", S, 3.4, 1.8), row("DNA-6", S, 2.9, 2.7), row("DNA-7", S, 3.5, 3.1),          # DNA: closest is DNA-6 (max 2.9)
                   row("Protein-5", S, 2.9, 2.7), row("Protein-8", S, 2.0, 2.0, fb=True),                       # Protein-8 is a weak pick: not used
                   row("Metabolite-5", S, 3.1, 3.5), row("RNA-1", "near SAM", 2.5, 4.6)])                     # RNA has no contact-both motif
got = [r["motif"] for r in choose(pl, S, closeness)]
assert got == ["DNA-1", "DNA-6", "Protein-5", "Metabolite-5"], got                       # panel order DNA, RNA, Protein, Metabolite; the RNA slot -> next overall closest (DNA-1: 3.4, before DNA-7: 3.5)
assert [r["why"] for r in choose(pl, S, closeness)][0] == "next overall"
assert "Protein-8" in [r["motif"] for r in choose(pl, S, closeness, allow_weak=True)]

D = "distal"
pd2 = pd.DataFrame([row("DNA-4", D, 5.4, 6.2), row("RNA-6", D, 4.6, 8.4), row("RNA-8", D, 9.1, 4.0), row("Metabolite-3", D, 16.4, 14.4), row("Protein-6", D, 12.5, 10.4, in_class=False),
                    row("Protein-5", S, 2.9, 2.7)])
got = [r["motif"] for r in choose(pd2, D, remoteness)]
assert got == ["DNA-4", "RNA-6", "RNA-8", "Metabolite-3"], got                            # RNA-6 (min 4.6) before RNA-8 (4.0); the weak Protein-6 is skipped, RNA-8 fills the slot
few = pd.DataFrame([row("DNA-4", D, 5.4, 6.2), row("RNA-1", "near SAM", 2.5, 4.6), row("Protein-5", S, 2.9, 2.7)])
got = choose(few, D, remoteness)
assert [r["motif"] for r in got][0] == "DNA-4" and len(got) == 3 and got[1]["why"].startswith("next overall, other placement")      # fewer than four: filled from the others, flagged
sel = select(pd.concat([pl, pd2[pd2.motif != "Protein-5"]], ignore_index=True))
assert sel.groupby("figure").size().to_dict() == {9: 4, 10: 4}
print("select_structures OK")
