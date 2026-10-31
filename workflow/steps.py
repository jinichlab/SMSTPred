"""Registry of every runnable step, in dependency order.

Each step runs from its own folder (`cwd`, relative to the repository root) with `cmd`. `make verify` runs the steps in a scratch
copy of the repository and compares every file a step regenerates with the committed version (see compare.py).

tier A   CPU, from the archived Zenodo data and this repository only (minutes)
tier B   GPU re-runs (`make verify-gpu`): outcomes are compared with stated tolerances by workflow/gpu_checks.py (field `check`)
tier C   needs a download or an external install (see REPRODUCING.md)

Fields: id, stage ("1".."7" or "figures"), tier, cwd, cmd, paper (what it reproduces), timeout (s), rtol (numeric tolerance),
png_tol (mean absolute pixel difference allowed, out of 255), compare_to ([(new, committed)] pairs when a step writes elsewhere),
check (name of a function in gpu_checks.py that replaces the file-by-file comparison).
"""


def S(id, stage, cwd, cmd, paper, tier="A", timeout=1800, rtol=1e-6, png_tol=12.0, compare_to=None, needs=(), check=None):
    return dict(id=id, stage=stage, tier=tier, cwd=cwd, cmd=cmd, paper=paper, timeout=timeout, rtol=rtol,
                png_tol=png_tol, compare_to=compare_to or [], needs=list(needs), check=check)


ORG = "Dickeya_undicola"          # one small organism, used to check the deployment model end to end

STEPS = [
    # ---------------------------------------------------------------- Stage 1: clustering
    S("s1_clustering", "1", "1_clustering", "python sam_clustering.py",
      "Stage 1 clustering tables at thresholds 0.5-0.95; Fig S1", timeout=7200, rtol=1e-5,
      needs=["prediction_files_mlp", "tree_of_life_ProtT5"]),
    # ---------------------------------------------------------------- Stage 2: benchmark statistics, tables, supplementary figures
    S("s2_bakeoff_stats", "2", "2_training", "python classifier_mcc_stats.py",
      "Friedman / Bonferroni-Dunn / Holm-Wilcoxon statistics of the benchmark (stats_report.txt)"),
    S("s2_fig2a_zoom", "2", "2_training", "python plot_cluster_mcc_zoom.py", "Fig 2A (cluster MCC, Friedman, Bonferroni-Dunn) and Fig S4"),
    S("s2_sam_vs_up_stats", "2", "2_training", "python sam_vs_up_prostt5_stats.py", "Fig 2B statistics (Nadeau-Bengio paired test)"),
    S("s2_sam_vs_up_plot", "2", "2_training", "python plot_sam_vs_up_prostt5.py", "Fig 2B bars"),
    S("s2_fig2_main", "2", "2_training", "python plot_fig2_main.py", "Fig 2 (A + B)"),
    S("s2_table1", "2", "2_training", "python build_metrics_table.py", "Table 1"),
    S("s2_figs_S2_S5", "2", "2_training", "python plot_mcc_curves_grid.py", "Figs S2 and S5 (MCC versus epoch)"),
    S("s2_figs_S3_S6", "2", "2_training", "python plot_train_val_loss_curves.py", "Figs S3 and S6 (train / validation loss)",
      needs=["train_val_loss_curves"]),
    S("s2_supp_stats", "2", "2_training", "python benchmark_supplement_stats.py",
      "Table S4 (Holm-adjusted Wilcoxon), Friedman tests for Fig 2A / Fig S4, the eleven ordinary iterations, per-iteration summary"),
    # ---------------------------------------------------------------- Stage 3: deployment-model inference
    S("s3_predict", "3", "3_predictions",
      f"python predict.py --input ../data/tree_of_life_ProstT5/{ORG}_embeddings.h5 --out-dir _verify_predict",
      "Deployment model: scores, calls, ensemble disagreement and applicability-domain flag for one organism",
      compare_to=[(f"3_predictions/_verify_predict/{ORG}_embeddings_predictions.csv",
                   f"3_predictions/predictions_output/{ORG}_embeddings_predictions.csv")],
      needs=["tree_of_life_ProstT5"], rtol=1e-4),
    S("s3_summary", "3", "3_predictions", "python summarize_predictions.py",
      "Prediction summary: 112,019 proteins, single / multi / no call, applicability domain, ensemble disagreement"),
    # ---------------------------------------------------------------- Stage 4: phylogenetic trends
    S("s4_lineage_stats", "4", "4_phylo_trends", "python lineage_statistics.py",
      "Lineage comparisons in the text (Micrococcales, Chlorophyta, plants, Actinomycetota; Mann-Whitney)"),
    # ---------------------------------------------------------------- Stage 5-6: UMAP, dark clusters, composite figures
    S("s5_fig6_umap", "5", "5_umap", "python plot_umap_static.py", "Fig 6 (static UMAP)"),
    S("s5_dark_definition", "5", "5_umap", "python dark_cluster_definition.py",
      "Dark-cluster definition: 868 clusters, 510 without a reviewed member, 33 dark (3,478 proteins, 3.6%)"),
    S("s6_dark_table", "6", "6_dark_cluster_annotation", "python build_table2.py", "Dark-cluster table (data behind Tables S1 and S2)"),
    S("s6_horizyn_direction", "6", "6_dark_cluster_annotation/horizyn", "python validation_by_direction.py",
      "Horizyn validation: top-1 is a SAM-methylation reaction, 44.0% forward only vs 85.7% either direction (chance 2.3% vs 4.7%)",
      needs=["horizyn_dark_scores", "horizyn_rhea_reaction_embeddings"]),
    # ---------------------------------------------------------------- Stage 7: motif analysis
    S("s7_domain_enrichment", "7", "7_motif_analysis", "python 3_domain_enrichment.py", "InterPro domain enrichment (Table S3 data)",
      needs=["motif_analysis_intermediates"]),
    S("s7_motif_table", "7", "7_motif_analysis", "python build_motif_table.py", "Table 2 (26 class-specific motifs)", needs=["motif_analysis_intermediates"]),
    S("s7_motif_status", "7", "7_motif_analysis", "python fimo_all_manifest.py", "Per-protein motif status (12,866 of 106,608 motif-less)",
      needs=["motif_analysis_intermediates"]),
    S("s7_motif_specificity", "7", "7_motif_analysis", "python fimo_cross_manifest.py", "Motif class specificity (26 of 32 at least 2-fold)",
      needs=["motif_analysis_intermediates"]),
    S("s7_family_reps", "7", "7_motif_analysis", "python 2_make_reps.py", "50%-identity family representatives (22,976)", timeout=3600,
      needs=["motif_analysis_intermediates"]),
    S("s7_fig11_nomotif", "7", "7_motif_analysis", "python no_motif_umap.py", "Fig 11 (motif-less enzymes on the UMAP)",
      needs=["motif_analysis_intermediates"]),
    # ---------------------------------------------------------------- figure and table builders (order matters)
    S("f_fig3_widen", "figures", "figures/Fig3_tree_landscape", "python widen_rings.py", "Fig 3, step 1 (widen the iTOL gradient rings)"),
    S("f_fig3_tree", "figures", "figures/Fig3_tree_landscape", "python make_tree_figure.py", "Fig 3, step 2 (tree + legend)"),
    S("f_fig1", "figures", "figures/Fig1_overview_schematic", "python build_fig1_flow.py", "Fig 1 (study overview)"),
    S("f_fig4_panels", "figures", "figures/Fig4_trends", "python build_fig4_panels.py", "Fig 4A-D"),
    S("f_fig4_combined", "figures", "figures/Fig4_trends", "python make_fig4_trends.py", "Fig 4 composite"),
    S("f_fig5", "figures", "figures/Fig5_proteome_scaling", "python proteome_vs_sam_domain_highlight.py", "Fig 5 (proteome scaling)"),
    S("f_fig7_grids", "figures", "figures/Fig7_clusters", "python build_cluster_grids.py", "Fig 7A-D cluster grids"),
    S("f_fig7_combined", "figures", "figures/Fig7_clusters", "python make_fig7_clusters.py", "Fig 7 composite"),
    S("f_fig8_reaction_profiles", "figures", "figures/Fig8_reaction_profiles", "python make_fig8_reaction_profiles.py",
      "Fig 8 (Horizyn reaction profiles of four dark clusters)", needs=["horizyn_dark_scores"]),
    S("f_fig9_structures", "figures", "figures/Fig9_10_structures", "python make_panels.py", "Structure panels (Figs 9 and 10)"),
    S("f_table_S1", "figures", "figures/TableS1_actinomycetota_dark_clusters", "python build_table.py", "Table S1"),
    S("f_table_S2", "figures", "figures/TableS2_plant_dark_clusters", "python build_table.py", "Table S2"),
    S("f_table_S3", "figures", "figures/TableS3_domain_enrichment", "python build_domain_enrichment_table.py", "Table S3"),

    # ---------------------------------------------------------------- tier B: GPU re-runs (the run_*.sh drivers, unchanged)
    S("b_deployment", "2", "2_training", "bash run_train_deployment.sh",
      "Deployment model re-trained from scratch: thresholds, out-of-fold MCC (0.969), out-of-fold scores", tier="B", check="deployment", timeout=3600),
    S("b_predict_all", "3", "3_predictions", "bash run_predict.sh",
      "Predictions for all 1,259 proteomes (112,019 proteins) from the re-trained model", tier="B", check="predictions", timeout=3600,
      needs=["tree_of_life_ProstT5"]),
    S("b_bakeoff_mlp", "2", "2_training", "bash run_bakeoff_mlp_only.sh",
      "Bake-off: the six MLP variants over the 16 cluster-disjoint iterations (Fig 2A, Table 1 rows)", tier="B", check="bakeoff", timeout=7200),
    S("b_loss_curves", "2", "2_training",
      "rm -rf train_val_loss_output && if command -v sbatch > /dev/null; then sbatch --wait run_train_val_loss.sh; else bash run_train_val_loss.sh; fi",
      "Train / validation loss curves behind Figs S3 and S6 (ten sweeps, SLURM array when sbatch exists)", tier="B", check="loss_curves", timeout=21600),
]

BY_ID = {s["id"]: s for s in STEPS}
