"""up_splits: clusters with a test protein are removed whole, validation is val_frac of the REMAINING clusters, every cluster is in exactly one partition, the split is deterministic and does not depend on how many are made."""
import numpy as np
import pandas as pd

import up_splits as u


def table(n_clusters=200, seed=0):
    rng = np.random.default_rng(seed)
    sizes = rng.integers(1, 6, n_clusters)
    rows = [(f"P{c}_{k}", c) for c, n in enumerate(sizes) for k in range(n)]
    return pd.DataFrame(rows, columns=["Entry", "cluster"]).assign(Metabolite=1, Protein=0, RNA=0, DNA=0)


def test_removed_validation_train_partition():
    up = table()
    test = ["P3_0", "P3_1", "P10_0", "P57_0"]                      # clusters 3, 10, 57
    tr, va, rm = u.make_up_split(up, test, 0.15, 42, 1)
    assert rm == {3, 10, 57}                                       # the WHOLE clusters of the test proteins, nothing else
    assert not (tr & va) and not (tr & rm) and not (va & rm) and (tr | va | rm) == set(up.cluster)
    assert len(va) == round(0.15 * (200 - 3))                      # 15 % of the REMAINING clusters


def test_deterministic_and_per_split_seed():
    up = table()
    a, b = u.make_up_split(up, ["P3_0"], 0.15, 42, 4), u.make_up_split(up, ["P3_0"], 0.15, 42, 4)
    assert a == b
    assert u.make_up_split(up, ["P3_0"], 0.15, 42, 5)[1] != a[1]   # another split index: another draw
    assert u.make_up_split(up, ["P3_0"], 0.15, 43, 4)[1] != a[1]   # another seed: another draw


def test_summary_counts():
    up = table(20)
    tr, va, rm = u.make_up_split(up, ["P1_0"], 0.25, 1, 1)
    s = u.summarize(up, [(tr, va, rm)]).set_index("partition")
    assert s.proteins.sum() == len(up) and s.loc["removed", "clusters"] == 1 and s.loc["train", "Metabolite"] == s.loc["train", "proteins"]


if __name__ == "__main__":
    for f in (test_removed_validation_train_partition, test_deterministic_and_per_split_seed, test_summary_counts):
        f()
        print("ok", f.__name__)
