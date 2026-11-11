#!/usr/bin/env python3
"""Variability of STREME across seeds: per variant the number of motifs with ratio_vs_bg >= 2 (the paper's 'class-discriminative' cut), by class, the consensus of every motif, and how many of a variant's motifs
reappear (same consensus up to 3 mismatches) in the reference run of the same inputs."""
import pandas as pd
from pathlib import Path
rows, cons = [], {}
for name, src, seed in (l.split() for l in open("variants.txt")):
    f = Path("runs") / name / "motif_specificity_matrix.csv"
    if not f.exists():
        rows.append(dict(variant=name, seed=seed, note="no matrix")); continue
    m = pd.read_csv(f)
    ok = m[m.ratio_vs_bg >= 2]
    cons[name] = m.set_index("motif").consensus.to_dict()
    r = dict(variant=name, seed=int(seed), n_motifs=len(m), n_ge2fold=len(ok))
    for c in ("DNA", "RNA", "Protein", "Metabolite"):
        r[c] = int((ok.own_class == c).sum())
    rows.append(r)
d = pd.DataFrame(rows); print(d.to_string(index=False))
def near(a, b):
    n = min(len(a), len(b))
    return abs(len(a) - len(b)) <= 2 and sum(x != y for x, y in zip(a[:n], b[:n])) <= 3
ref = {"arch": "arch_seed0", "rebuilt": "rebuilt_seed0"}
out = []
for name, c in cons.items():
    base = ref[name.split("_")[0]]
    if name == base or base not in cons: continue
    rc = list(cons[base].values())
    out.append(dict(variant=name, vs=base, motifs_found_in_reference=sum(any(near(v, x) for x in rc) for v in c.values()), of=len(c)))
o = pd.DataFrame(out); print(o.to_string(index=False))
d.to_csv("variability_counts.csv", index=False); o.to_csv("variability_overlap.csv", index=False)
pd.DataFrame([dict(variant=n, motif=k, consensus=v) for n, c in cons.items() for k, v in c.items()]).to_csv("variability_motifs.csv", index=False)
