#!/usr/bin/env python3
"""motif_structure_pick: placement classes, locating a motif in a chain (with one mismatch), and the pass/fail gates on a tiny synthetic structure."""
import tempfile
from pathlib import Path

import gemmi
from motif_structure_pick import evaluate, locate_all, majority_placement, placement, rank_key

assert placement(3.0, 2.0) == "near substrate/SAM" and placement(3.0, 9.0) == "near SAM" and placement(9.0, 3.9) == "near substrate" and placement(4.1, 4.1) == "distal"


def pdb_text(sam_x, lig_x):
    """Chain A = GGKSAV (6 residues, CA 4 A apart along x), SAH on chain A at x = sam_x, ligand XYZ at x = lig_x."""
    L, n = [], 1
    for i, rn in enumerate(["GLY", "GLY", "LYS", "SER", "ALA", "VAL"]):
        L.append(f"ATOM  {n:5d}  CA  {rn} A{i + 1:4d}    {4.0 * i:8.3f}{0:8.3f}{0:8.3f}  1.00 20.00           C")
        n += 1
    for name, x, res in (("N9", sam_x, "SAH"), ("C1", lig_x, "XYZ")):
        L.append(f"HETATM{n:5d}  {name:<3s} {res} A{100 + n:4d}    {x:8.3f}{3:8.3f}{0:8.3f}  1.00 20.00           C")
        n += 1
    return "\n".join(L) + "\nEND\n"


def run(sam_x, lig_x, matched="GGKSAV"):
    f = Path(tempfile.mkdtemp()) / "t.pdb"
    f.write_text(pdb_text(sam_x, lig_x))
    return evaluate(f, matched, "small_molecule", {"XYZ"})


st = gemmi.read_structure(str(Path(tempfile.mkdtemp()) / "x.pdb")) if False else None
f = Path(tempfile.mkdtemp()) / "t.pdb"; f.write_text(pdb_text(5.0, 6.0))
model = gemmi.read_structure(str(f))[0]
assert locate_all(model, "GGKSAV")[0][1:3] == (1, 6)
assert locate_all(model, "GGKSAA")[0][1:3] == (1, 6)           # one mismatch allowed
assert not locate_all(model, "GGWWAA")                          # two mismatches: not found

ev, why = run(5.0, 6.0)                                         # SAM and the ligand next to the motif and to each other
assert ev is not None and why == "ok" and ev["d_ss"] < 6.0, (ev, why)
ev, why = run(5.0, 30.0)                                        # ligand far from the SAM: not the active site (and far from the motif)
assert ev is None, why
ev, why = run(5.0, 6.0, "WWWWWW")
assert ev is None and "motif not found" in why
# ranking of several passing structures of one motif
def P(**k):
    d = dict(kind_fallback=False, in_class=True, placement="near SAM", reviewed=True, cofactor="SAH", resolution=2.0, fimo_p=1e-6)
    d.update(k)
    return d
typical, odd = P(), P(placement="distal")
maj = majority_placement([typical, P(), odd])
assert maj == "near SAM"
order = lambda *ps: sorted(ps, key=lambda q: rank_key(q, maj))
assert order(odd, typical)[0] is typical                                            # the placement most structures show beats an unusual one
assert order(P(in_class=False), P())[0]["in_class"]                                   # own class first
assert order(P(kind_fallback=True, in_class=True), P(in_class=False))[0]["in_class"] is False       # an own-kind substrate beats an own-class carrier with a fallback substrate
assert order(P(reviewed=False), P())[0]["reviewed"]
assert order(P(cofactor="SAH"), P(cofactor="SAM"))[0]["cofactor"] == "SAM"           # intact SAM before the product
assert order(P(resolution=2.9), P(resolution=1.6))[0]["resolution"] == 1.6
assert order(P(resolution=2.0, fimo_p=1e-3), P(resolution=2.2, fimo_p=1e-9))[0]["fimo_p"] == 1e-9   # same 0.5 A bucket... 2.0 vs 2.2 round to 2.0 and 2.0: the better match wins
assert order(P(resolution=None), P(resolution=3.5))[0]["resolution"] == 3.5          # unknown resolution last
assert majority_placement([P(placement="distal", in_class=False), P(placement="near SAM")]) == "near SAM"      # out-of-class ones do not vote when own-class ones exist
print("motif_structure_pick OK")
