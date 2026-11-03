#!/usr/bin/env python3
"""select_panels: the Fig 8 rule on a small made-up release, plus the real mapping of Rhea entries to their master reaction."""
import pandas as pd
from select_panels import family_ecs, master_map, select

# master 100 has entries 101 (LR), 102 (RL), 103 (BI); master 205 (not a multiple of 4!) has 206, 207, 208
directions = pd.DataFrame({"RHEA_ID_MASTER": [100, 205], "RHEA_ID_LR": [101, 206], "RHEA_ID_RL": [102, 207], "RHEA_ID_BI": [103, 208]})
m = master_map(directions)
assert m[102] == 100 and m[207] == 205 and m[205] == 205

rhea2ec = pd.DataFrame({"RHEA_ID": ["100", "205"], "DIRECTION": ["UN", "UN"], "MASTER_ID": ["100", "205"], "ID": ["2.1.1.4", "2.1.1.22"]})
rxn = pd.DataFrame({"rhea_lr_id": [101, 102, 206, 999], "reaction_cluster": [7, 7, 8, 8]})     # 999 has no EC; family 8 holds 2.1.1.22
fe = family_ecs(rxn, rhea2ec, m)
assert fe == {7: {"2.1.1.4"}, 8: {"2.1.1.22"}}, fe

clean = pd.DataFrame({"cluster": [1, 2, 3, 4, 5], "n": [50, 60, 70, 80, 90],
                      "top_ec": ["EC:2.1.1.4", "EC:2.1.1.22", "EC:2.1.1.4", "EC:2.1.1.4", "EC:2.1.1.22"],
                      "confidence": ["reliable (d=3.0)", "reliable (d=5.5)", "marginal (d=7.0)", "reliable (d=4.0)", "reliable (d=2.0)"],
                      "mtase": ["yes", "yes", "yes", "No", "yes"]})
fam = pd.DataFrame({"cluster": [1, 2, 3, 4, 5], "coherent_dominant_family": [7, 7, 7, 7, 8]})
# 1: reliable MTase, EC in family -> shown. 2: EC 2.1.1.22 not in family 7 -> no. 3: marginal -> no. 4: not an MTase -> no. 5: shown.
assert select(clean, fam, rxn, rhea2ec, m) == [1, 5]
assert select(clean, fam, rxn, rhea2ec, m, max_panels=1) == [5]      # more than max_panels: the largest
print("select_panels OK")
