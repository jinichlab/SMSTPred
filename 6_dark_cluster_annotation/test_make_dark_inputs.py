#!/usr/bin/env python3
"""make_dark_inputs.build: members of the dark clusters only, ordered as dark_clusters.csv, lineage attached; a missing sequence stops the run."""
import pandas as pd
from make_dark_inputs import build

dark = pd.DataFrame({"cluster": [9, 3], "lineage": ["Plant", "Actinomycetota"]})
coords = pd.DataFrame({"Entry": ["a", "b", "c", "d"], "Organism": ["x", "y", "z", "w"], "cluster": [3, 9, 5, 9]})
seqs = {"a": "MK", "b": "MA", "c": "MG", "d": "ML"}
m = build(dark, coords, seqs)
assert m["Entry"].tolist() == ["b", "d", "a"], m                 # cluster 9 first (as in dark_clusters.csv), c (cluster 5) is not dark
assert m["lineage"].tolist() == ["Plant", "Plant", "Actinomycetota"]
try:
    build(dark, coords, {"a": "MK"})
    raise AssertionError("a dark protein without a sequence must stop the run")
except SystemExit:
    pass
print("make_dark_inputs OK")
