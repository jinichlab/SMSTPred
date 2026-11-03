"""Names for Rhea reactions, from the equation text (<release>/reference_tables/rhea_equations.tsv, 0_uniprot_pull/fetch_reference_tables.py).

reaction_title("N-acetylserotonin + S-adenosyl-L-methionine = melatonin + S-adenosyl-L-homocysteine + H(+)") -> "N-acetylserotonin → melatonin": the cofactors and protons dropped,
the side holding SAM first. Used for the Fig 8 panel titles and the descriptions of the dark-cluster tables.
"""
import re

import pandas as pd

COFACTORS = {"S-adenosyl-L-methionine", "S-adenosyl-L-homocysteine", "H(+)", "H2O"}
SPLIT = re.compile(r"\s(?:<=>|=>|=|<=)\s")


def load_equations(path):
    """{master or entry id (int): equation text}"""
    d = pd.read_csv(path, sep="\t", dtype={"RHEA_ID": int}, quoting=3)
    return dict(zip(d["RHEA_ID"], d["equation"]))


def _side(s):
    out = []
    for t in s.split(" + "):
        t = re.sub(r"^\d+\s+", "", t.strip())
        if t and t not in COFACTORS:
            out.append(t)
    return out


def reaction_title(eq, max_len=60):
    """'substrate → product' of an equation (SAM side first); None if nothing is left once the cofactors are removed."""
    parts = SPLIT.split(eq)
    if len(parts) != 2:
        return None
    left, right = parts
    if "S-adenosyl-L-methionine" in right.split(" + ") and "S-adenosyl-L-methionine" not in left.split(" + "):
        left, right = right, left
    a, b = _side(left), _side(right)
    if not a or not b:
        return None
    t = f"{' + '.join(a)} → {' + '.join(b)}"
    return t if len(t) <= max_len else t[:max_len - 1] + "…"
