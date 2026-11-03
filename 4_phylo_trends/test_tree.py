#!/usr/bin/env python3
"""Tests of the release-mode tree: combine_predictions, reference_tree.collect_release and the in-house drawing (synthetic data, no release needed).  python test_tree.py"""
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / "3_predictions"), str(HERE.parent / "figures" / "Fig3_tree_landscape")]
import combine_predictions as cp
import draw_tree as dt
import reference_tree as rt


def fake_release(tmp):
    (tmp / "predictions").mkdir()
    (tmp / "embedding_inputs").mkdir()
    cols = "Entry,score_Metabolite,score_Protein,score_RNA,score_DNA,call_Metabolite,call_Protein,call_RNA,call_DNA,disagreement,ad_maxsim,in_domain\n"
    # proteome A: 4 proteins (2 RNA, 1 DNA, 1 none); proteome B: 2 proteins (1 Metabolite, 1 Protein+RNA)
    (tmp / "predictions" / "A_embeddings_predictions.csv").write_text(cols + "a1,0,0,1,0,0,0,1,0,0,0.9,1\na2,0,0,1,0,0,0,1,0,0,0.9,0\na3,0,0,0,1,0,0,0,1,0,0.9,1\na4,0,0,0,0,0,0,0,0,0,0.5,0\n")
    (tmp / "predictions" / "B_embeddings_predictions.csv").write_text(cols + "b1,1,0,0,0,1,0,0,0,0,0.9,1\nb2,0,1,1,0,0,1,1,0,0,0.9,1\n")
    pd.DataFrame({"proteome": ["A"] * 4 + ["B"] * 2, "Entry": ["a1", "a2", "a3", "a4", "b1", "b2"], "Length": [100, 200, 300, 400, 500, 600]}).to_csv(tmp / "embedding_inputs" / "proteome_members.csv", index=False)
    pd.DataFrame({"proteome_id": ["UP1", "UP2", "UP3"], "stem": ["A", "B", "C"], "status": ["ok", "ok", "empty"], "n_all": [1000, 500, 0], "n_ipr_all": [5, 3, 0]}).to_csv(tmp / "proteome_summary.csv", index=False)
    return tmp


def test_combine_and_collect():
    with tempfile.TemporaryDirectory() as t:
        rel = fake_release(Path(t))
        out = cp.combine(rel)
        assert len(out) == 6 and out.set_index("Entry").Length["b2"] == 600 and set(out.proteome) == {"A", "B"}
        out.to_csv(rel / "predictions_all_proteomes.csv.gz", index=False)
        df = pd.DataFrame({"Proteome Id": ["UP1", "UP2", "UP3"], "tree_id": ["org_A", "org_B", "org_C"]})
        stats = defaultdict(lambda: {"total": 0, "sams": 0, "predicted": 0, "RNA": 0, "DNA": 0, "Protein": 0, "Metabolite": 0})
        rt.collect_release(df, rel, stats)
        a, b = stats["org_A"], stats["org_B"]
        assert (a["total"], a["sams"], a["predicted"], a["RNA"], a["DNA"], a["Protein"], a["Metabolite"]) == (1000, 5, 4, 2, 1, 0, 0)
        assert (b["sams"], b["predicted"], b["RNA"], b["Protein"], b["Metabolite"]) == (3, 2, 1, 1, 1), "a multi-label protein counts once in each label it is called"
        assert "org_C" not in stats, "an empty proteome is not in the tree"


def test_newick_layout_and_draw():
    root = dt.parse_newick("((a:1,b:1)X:1,(c:1,d:1,e:1)Y:1)Z:1;")
    leaves = dt.layout(root, 0.5)
    assert [l.name for l in leaves] == list("abcde")
    angles = [l.angle for l in leaves]
    assert all(x > y for x, y in zip(angles, angles[1:])), "clockwise from the top"
    assert abs(angles[0] - (3.14159265 / 2 - 2 * 3.14159265 / 10)) < 1e-6 and all(l.radius == 1.0 for l in leaves)
    with tempfile.TemporaryDirectory() as t:
        src = Path(t)
        (src / "reference_tree.nwk").write_text("((a:1,b:1)X:1,(c:1,d:1,e:1)Y:1);")
        for stem in ("itol_pct_protein", "itol_pct_dna", "itol_pct_rna", "itol_pct_metabolite", "itol_pct_sams"):
            (src / f"{stem}.txt").write_text("DATASET_GRADIENT\nSEPARATOR TAB\nDATASET_LABEL\tx\nCOLOR\t#238b45\nCOLOR_MIN\t#f7fbff\nCOLOR_MAX\t#238b45\nDATA\n" + "".join(f"{k}\t{v}\n" for k, v in zip("abcde", (0.0, 0.1, 0.2, 0.3, 0.4))))
        (src / "itol_n_sams.txt").write_text("DATASET_SIMPLEBAR\nSEPARATOR TAB\nDATASET_LABEL\tn\nCOLOR\t#2b8cbe\nDATA\n" + "".join(f"{k}\t{v}\n" for k, v in zip("abcde", (20, 30, 40, 50, 60))))
        (src / "itol_clade_highlights.txt").write_text("TREE_COLORS\nSEPARATOR TAB\nDATA\na|b\trange\t#33a02c\tPlants\nc|e\trange\t#1f78b4\tMetazoa\n")
        dt.build(src, src, "fig", 0.5)
        assert (src / "fig.png").stat().st_size > 5000 and (src / "fig.svg").exists() and (src / "fig.pdf").exists()


if __name__ == "__main__":
    test_combine_and_collect()
    test_newick_layout_and_draw()
    print("all tree tests passed")
