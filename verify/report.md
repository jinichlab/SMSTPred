Last verification: 2026-11-11T18:20:06 (Python 3.12.15, numpy 2.1.3, pandas 2.2.3, scipy 1.15.2, torch 2.8.0, matplotlib 3.10.9). 14 of 18 steps reproduce a committed result, 4 build their output without a committed reference to compare, 0 fail.

| Step | Reproduces | Result | Files compared | Time |
|---|---|---|---|---|
| `t_binning` | Substrate typing rules (Rhea compound classes -> Metabolite / Protein / RNA / DNA) | built (no reference) | ran; wrote no files | 1 s |
| `t_training` | Cluster-disjoint splits, bake-off, MTase-vs-UP and ensemble code, MCC threshold search | built (no reference) | ran; wrote no files | 61 s |
| `t_annotation` | Cluster selection, dark-protein inputs, CLEAN summary, Horizyn family table, Fig 8 rule, motif structure search and pick | built (no reference) | ran; wrote no files | 5 s |
| `s1_cluster` | 502 clusters of the benchmark set, 16 splits, cluster columns of the label table | built (no reference) | 2 new-file | 19 s |
| `s2_bakeoff_report` | Fig 2A statistics, Table 1, the paired tests against the deployed model, the train / validation loss figures | PASS | 18 identical, 1 image-match | 15 s |
| `s2_sam_vs_up_report` | Fig 2B: MTase-trained vs UP-trained, Wilcoxon-Holm and Nadeau-Bengio | PASS | 8 identical | 4 s |
| `s2_fig2` | Fig 2 composite (A + B) | PASS | 2 identical | 6 s |
| `s3_predict` | Scores, calls, ensemble disagreement and applicability-domain flag for every proteome protein | PASS | 1249 new-file, 1 numeric-match | 28 s |
| `s4_tree` | Fig 3: taxonomy tree and per-organism annotation | PASS | 10 identical, 1 numeric-match, 2 skipped | 10 s |
| `s4_fig4` | Fig 4 and its Mann-Whitney / Kruskal-Wallis statistics | PASS | 7 identical | 4 s |
| `s4_fig5` | Fig 5 power-law fits, level differences, clade-bootstrap CIs, lineage statistics | PASS | 11 identical, 1 image-match, 2 skipped | 35 s |
| `s5_umap` | The map: the archived UMAP coordinates and clusters are reused (frozen), scores and calls refreshed from the predictions; the static map is redrawn | PASS | 2 identical, 1 numeric-match, 1 skipped | 42 s |
| `s5_clusters` | Dark / metazoan-protein / vertebrate-DNA cluster selection and Fig 7 | PASS | 11 identical | 5 s |
| `s6_tables` | Horizyn family table, CLEAN summary, Tables S1 and S2, Fig 7 with the Fig 8 clusters | PASS | 17 identical | 16 s |
| `s6_fig8` | Fig 8: Horizyn reaction profiles of the rule-chosen dark clusters | PASS | 2 identical | 6 s |
| `s7_no_motif_umap` | Fig 11: enzymes with no discriminative motif on the map | PASS | 2 identical | 2 s |
| `s7_motif_pick` | Motif structure pick by rule, Table 2 | PASS | 77 identical, 185 new-file | 64 s |
| `s9_key_numbers` | key_numbers.md (every number quoted in the text) and the figure index | PASS | 1 identical, 24 new-file | 3 s |

Last verification: 2026-11-11T18:20:06 (Python 3.12.15, numpy 2.1.3, pandas 2.2.3, scipy 1.15.2, torch 2.8.0, matplotlib 3.10.9). 1 of 1 steps reproduce a committed result, 0 build their output without a committed reference to compare, 0 fail.

| Step | Reproduces | Result | Checks | Time |
|---|---|---|---|---|
| `b_deploy` | Deployment ensemble re-trained from the labels: thresholds, out-of-fold MCC, out-of-fold scores | PASS | 8 of 8 passed | 93 s |
| &nbsp;&nbsp;`b_deploy` check | same proteins, clusters and folds | PASS | 7582 proteins | |
| &nbsp;&nbsp;`b_deploy` check | six ensemble members written | PASS |  | |
| &nbsp;&nbsp;`b_deploy` check | every label called at the same fixed threshold as archived | PASS | {'Metabolite': 0.5, 'Protein': 0.5, 'RNA': 0.5, 'DNA': 0.5} | |
| &nbsp;&nbsp;`b_deploy` check | overall out-of-fold MCC within 0.01 | PASS | 0.9549 vs 0.9551 | |
| &nbsp;&nbsp;`b_deploy` check | per-label out-of-fold MCC within 0.02 | PASS | Metabolite 0.962 vs 0.962, Protein 0.913 vs 0.913, RNA 0.968 vs 0.969, DNA 0.976 vs 0.976 | |
| &nbsp;&nbsp;`b_deploy` check | out-of-fold scores within 0.1 (max over proteins) | PASS | max difference 0.0094 | |
| &nbsp;&nbsp;`b_deploy` check | mean out-of-fold score difference within 0.01 | PASS | mean difference 0.0001 | |
| &nbsp;&nbsp;`b_deploy` check | applicability-domain cutoff within 0.01 | PASS | 0.8130 vs 0.8130 | |
