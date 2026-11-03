"""tvt_splits: every cluster in exactly one partition, the sizes follow the fractions, same seed = same splits, the first k splits do not depend on how many are made."""
import numpy as np
import pandas as pd

import tvt_splits as t


def test_partition_properties():
    ids = list(range(100))
    sp = t.make_tvt_splits(ids, 5, 0.2, 0.1, 7)
    for tr, va, te in sp:
        assert len(te) == 20 and len(va) == 10 and len(tr) == 70
        assert not (tr & va) and not (tr & te) and not (va & te) and (tr | va | te) == set(ids)
    assert len({frozenset(s[2]) for s in sp}) == 5                                # the splits differ from each other


def test_deterministic_and_prefix_stable():
    a, b = t.make_tvt_splits(range(50), 4, 0.2, 0.1, 42), t.make_tvt_splits(range(50), 4, 0.2, 0.1, 42)
    assert a == b
    assert t.make_tvt_splits(range(50), 2, 0.2, 0.1, 42) == a[:2]                 # 2 splits = the first 2 of 4
    assert t.make_tvt_splits(range(50), 2, 0.2, 0.1, 43) != a[:2]


def test_summary_counts():
    df = pd.DataFrame({"Entry": list("abcdef"), "cluster": [0, 0, 1, 2, 3, 3], "Metabolite": [1, 0, 0, 0, 1, 0], "Protein": [0, 1, 0, 0, 0, 0],
                       "RNA": [0, 0, 1, 1, 0, 1], "DNA": [0, 0, 0, 0, 0, 0]})
    s = t.summarize(df, [({0, 1}, {2}, {3})])
    r = s.set_index("partition")
    assert r.loc["train", "proteins"] == 3 and r.loc["train", "Metabolite"] == 1 and r.loc["train", "Protein"] == 1 and r.loc["train", "RNA"] == 1
    assert r.loc["validation", "proteins"] == 1 and r.loc["validation", "RNA"] == 1 and r.loc["test", "proteins"] == 2 and r.loc["test", "Metabolite"] == 1


if __name__ == "__main__":
    for f in (test_partition_properties, test_deterministic_and_prefix_stable, test_summary_counts):
        f()
        print("ok", f.__name__)
