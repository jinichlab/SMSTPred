Last verification: 2026-10-30T17:31:28 (Python 3.12.15, numpy 2.1.3, pandas 2.2.3, scipy 1.15.2, torch 2.8.0, matplotlib 3.10.9). 26 of 36 steps reproduce a committed result, 10 build their output without a committed reference to compare, 0 fail.

| Step | Reproduces | Result | Files compared | Time |
|---|---|---|---|---|
| `s1_clustering` | Stage 1 clustering tables at thresholds 0.5-0.95; Fig S1 | PASS | 2 identical, 2 image-match, 1 new-file, 2 numeric-match | 59 s |
| `s2_bakeoff_stats` | Friedman / Bonferroni-Dunn / Holm-Wilcoxon statistics of the benchmark (stats_report.txt) | PASS | 3 identical | 3 s |
| `s2_fig2a_zoom` | Fig 2A (cluster MCC, Friedman, Bonferroni-Dunn) and Fig S4 | PASS | 2 identical | 1 s |
| `s2_sam_vs_up_stats` | Fig 2B statistics (Nadeau-Bengio paired test) | PASS | 2 identical | 1 s |
| `s2_sam_vs_up_plot` | Fig 2B bars | PASS | 1 identical | 1 s |
| `s2_fig2_main` | Fig 2 (A + B) | PASS | 1 identical | 1 s |
| `s2_table1` | Table 1 | PASS | 1 identical | 1 s |
| `s2_figs_S2_S5` | Figs S2 and S5 (MCC versus epoch) | PASS | 2 identical | 3 s |
| `s2_figs_S3_S6` | Figs S3 and S6 (train / validation loss) | PASS | 12 identical | 31 s |
| `s2_supp_stats` | Table S4 (Holm-adjusted Wilcoxon), Friedman tests for Fig 2A / Fig S4, the eleven ordinary iterations, per-iteration summary | PASS | 6 identical | 1 s |
| `s3_predict` | Deployment model: scores, calls, ensemble disagreement and applicability-domain flag for one organism | PASS | 1 identical, 1 new-file | 4 s |
| `s3_summary` | Prediction summary: 112,019 proteins, single / multi / no call, applicability domain, ensemble disagreement | PASS | 2 identical | 3 s |
| `s4_lineage_stats` | Lineage comparisons in the text (Micrococcales, Chlorophyta, plants, Actinomycetota; Mann-Whitney) | PASS | 2 identical | 1 s |
| `s5_fig6_umap` | Fig 6 (static UMAP) | PASS | 1 identical | 2 s |
| `s5_dark_definition` | Dark-cluster definition: 868 clusters, 510 without a reviewed member, 33 dark (3,478 proteins, 3.6%) | PASS | 2 identical | 1 s |
| `s6_dark_table` | Dark-cluster table (data behind Tables S1 and S2) | PASS | 1 identical | 1 s |
| `s6_horizyn_direction` | Horizyn validation: top-1 is a SAM-methylation reaction, 44.0% forward only vs 85.7% either direction (chance 2.3% vs 4.7%) | PASS | 2 identical | 7 s |
| `s7_domain_enrichment` | InterPro domain enrichment (Table S3 data) | PASS | 1 identical | 2 s |
| `s7_motif_table` | Table 2 (26 class-specific motifs) | PASS | 2 identical | 2 s |
| `s7_motif_status` | Per-protein motif status (12,866 of 106,608 motif-less) | PASS | 2 identical | 1 s |
| `s7_motif_specificity` | Motif class specificity (26 of 32 at least 2-fold) | PASS | 2 identical, 1 new-file | 4 s |
| `s7_family_reps` | 50%-identity family representatives (22,976) | PASS | 5 identical | 1 s |
| `s7_fig11_nomotif` | Fig 11 (motif-less enzymes on the UMAP) | PASS | 2 identical | 1 s |
| `f_fig3_widen` | Fig 3, step 1 (widen the iTOL gradient rings) | built (no reference) | 1 new-file | 1 s |
| `f_fig3_tree` | Fig 3, step 2 (tree + legend) | built (no reference) | 2 new-file | 9 s |
| `f_fig1` | Fig 1 (study overview) | built (no reference) | 3 new-file | 2 s |
| `f_fig4_panels` | Fig 4A-D | built (no reference) | 4 new-file | 3 s |
| `f_fig4_combined` | Fig 4 composite | built (no reference) | 1 new-file | 1 s |
| `f_fig5` | Fig 5 (proteome scaling) | built (no reference) | 1 new-file | 1 s |
| `f_fig7_grids` | Fig 7A-D cluster grids | built (no reference) | 4 new-file | 2 s |
| `f_fig7_combined` | Fig 7 composite | built (no reference) | 1 new-file | 1 s |
| `f_fig8_reaction_profiles` | Fig 8 (Horizyn reaction profiles of four dark clusters) | built (no reference) | 1 new-file | 4 s |
| `f_fig9_structures` | Structure panels (Figs 9 and 10) | built (no reference) | 2 new-file | 1 s |
| `f_table_S1` | Table S1 | PASS | 1 identical, 1 new-file | 1 s |
| `f_table_S2` | Table S2 | PASS | 1 identical, 1 new-file | 1 s |
| `f_table_S3` | Table S3 | PASS | 1 identical, 1 new-file | 1 s |

Last verification: 2026-10-30T17:31:28 (Python 3.12.15, numpy 2.1.3, pandas 2.2.3, scipy 1.15.2, torch 2.8.0, matplotlib 3.10.9). 4 of 4 steps reproduce a committed result, 0 build their output without a committed reference to compare, 0 fail.

| Step | Reproduces | Result | Checks | Time |
|---|---|---|---|---|
| `b_deployment` | Deployment model re-trained from scratch: thresholds, out-of-fold MCC (0.969), out-of-fold scores | PASS | 4 of 4 passed | 13 s |
| `b_predict_all` | Predictions for all 1,259 proteomes (112,019 proteins) from the re-trained model | PASS | 4 of 4 passed | 32 s |
| `b_bakeoff_mlp` | Bake-off: the six MLP variants over the 16 cluster-disjoint iterations (Fig 2A, Table 1 rows) | PASS | 3 of 3 passed | 168 s |
| `b_loss_curves` | Train / validation loss curves behind Figs S3 and S6 (ten sweeps, SLURM array when sbatch exists) | PASS | 4 of 4 passed | 4206 s |
| &nbsp;&nbsp;`b_deployment` check | same proteins, clusters and folds | PASS | 7405 proteins | |
| &nbsp;&nbsp;`b_deployment` check | per-label thresholds within 0.02 | PASS | Metabolite 0.626 vs 0.634, Protein 0.635 vs 0.632, RNA 0.827 vs 0.826, DNA 0.109 vs 0.109 | |
| &nbsp;&nbsp;`b_deployment` check | per-label out-of-fold MCC within 0.01 | PASS | Metabolite 0.970 vs 0.970, Protein 0.971 vs 0.973, RNA 0.990 vs 0.990, DNA 0.945 vs 0.945 | |
| &nbsp;&nbsp;`b_deployment` check | out-of-fold scores within 0.03 | PASS | max difference 0.0117 | |
| &nbsp;&nbsp;`b_predict_all` check | same organism files | PASS | 1259 files | |
| &nbsp;&nbsp;`b_predict_all` check | identical class calls for at least 99.5% of proteins | PASS | 99.80% of 112,019 proteins | |
| &nbsp;&nbsp;`b_predict_all` check | scores within 0.05 | PASS | max difference 0.0228 | |
| &nbsp;&nbsp;`b_predict_all` check | applicability-domain flags identical | PASS |  | |
| &nbsp;&nbsp;`b_bakeoff_mlp` check | every MLP row of the re-run has a committed counterpart | PASS | 192 of 192 rows | |
| &nbsp;&nbsp;`b_bakeoff_mlp` check | max absolute MCC difference per iteration at most 0.06 | PASS | max 0.040, mean 0.0062 | |
| &nbsp;&nbsp;`b_bakeoff_mlp` check | mean MCC per method within 0.01 of the committed run | PASS | MLP_ESM-IF 0.9137 vs 0.9118; MLP_ProstT5 0.9563 vs 0.9514; MLP_ProstT5+ESM-IF 0.9483 vs 0.9474; MLP_ProtT5 0.9504 vs 0.9489; MLP_ProtT5+ESM-IF 0.9498 vs 0.9520; MLP_SaProt 0.9512 vs 0.9503 | |
| &nbsp;&nbsp;`b_loss_curves` check | a curves file for every committed sweep | PASS | 10 of 10 sweeps | |
| &nbsp;&nbsp;`b_loss_curves` check | median absolute validation-loss difference over all epochs at most 0.01 | PASS | median 0.0015, 95th percentile 0.1863 | |
| &nbsp;&nbsp;`b_loss_curves` check | median absolute final-epoch validation-loss difference at most 0.01 | PASS | median 0.0072, max 22.6662 | |
| &nbsp;&nbsp;`b_loss_curves` check | median absolute checkpoint-MCC difference at most 0.02 | PASS | median 0.0000, 95th percentile 0.0756 | |
