"""Registry of the steps `make verify` runs, in dependency order, for the release pipeline.

`make data` assembles the archived Zenodo record into a release folder (workflow/fetch_data.py: data/uniprot_pull/release_<tag>/). `make verify` copies that folder
to a scratch folder, runs each stage's own driver script on the copy (RELEASE_DIR points at it) and compares every file the stage wrote with the archived file of the
same name (compare.py). Steps whose results are not in the record are reported as BUILT (they ran and produced output, there is nothing to compare). The old registry of the
2024-era paper (iTOL tree, 0.6 embedding clustering, ProstT5 deployment) is kept as workflow/legacy_steps_paper_2024.py and is not run.

tier A   CPU, from the archived record and this repository only
tier B   GPU: re-train and compare outcomes with stated tolerances (workflow/gpu_checks.py, field `check`)
tier C   needs a download or an external install (UniProt pull, CLEAN, Horizyn, ChimeraX, ColabFold): documented in REPRODUCING.md, never part of `make verify`

Fields: id, stage, tier, cmd (run with bash from the scratch copy of the repository), paper (what it reproduces), timeout (s), rtol, png_tol, check, skip_compare (files
a step rewrites that must not be compared, e.g. logs).
"""


def S(id, stage, cmd, paper, tier="A", timeout=3600, rtol=1e-5, png_tol=12.0, check=None, skip_compare=()):
    return dict(id=id, stage=stage, tier=tier, cwd=".", cmd=cmd, paper=paper, timeout=timeout, rtol=rtol, png_tol=png_tol, compare_to=[], needs=[],
                check=check, skip_compare=list(skip_compare))


STEPS = [
    # ---------------------------------------------------------------- offline unit tests of the rules and statistics
    S("t_binning", "0", "cd 0_uniprot_pull && python test_substrate_binning.py", "Substrate typing rules (Rhea compound classes -> Metabolite / Protein / RNA / DNA)", timeout=900),
    S("t_training", "1", "cd 1_clustering && python test_tvt_splits.py && python test_up_splits.py && cd ../2_training && python test_bakeoff.py && python test_sam_vs_up.py && python test_deploy_ensemble.py && python test_threshold_search.py",
      "Cluster-disjoint splits, bake-off, MTase-vs-UP and ensemble code, MCC threshold search", timeout=1800),
    S("t_annotation", "5", "make annotation-test", "Cluster selection, dark-protein inputs, CLEAN summary, Horizyn family table, Fig 8 rule, motif structure search and pick", timeout=1800),
    # ---------------------------------------------------------------- Stage 1: MMseqs2 30% clustering and the 16 cluster-disjoint splits, re-derived from the archived labels
    S("s1_cluster", "1", "CLUSTER_WHAT=sam bash 1_clustering/run_mmseqs_clusters.sh", "502 clusters of the benchmark set, 16 splits, cluster columns of the label table", timeout=3600),
    # ---------------------------------------------------------------- Stage 2: benchmark statistics and figures from the archived per-split results
    S("s2_bakeoff_report", "2", "bash 2_training/run_bakeoff_report.sh", "Fig 2A statistics, Table 1, Table S4, the train / validation loss figures", timeout=1800),
    S("s2_sam_vs_up_report", "2", "bash 2_training/run_sam_vs_up_report.sh", "Fig 2B: MTase-trained vs UP-trained, Wilcoxon-Holm and Nadeau-Bengio", timeout=1800),
    S("s2_fig2", "2", "bash 2_training/run_fig2.sh", "Fig 2 composite (A + B)", timeout=900),
    # ---------------------------------------------------------------- Stage 3: the deployed ensemble applied to every reference proteome (CPU inference)
    S("s3_predict", "3", "bash 3_predictions/run_predict_release.sh", "Scores, calls, ensemble disagreement and applicability-domain flag for every proteome protein", timeout=7200, rtol=1e-3),
    # ---------------------------------------------------------------- Stage 4-5: tree, trends, scaling, UMAP clusters
    S("s4_tree", "4", "bash 4_phylo_trends/run_tree.sh", "Fig 3: taxonomy tree and per-organism annotation", timeout=1800),
    S("s4_fig4", "4", "bash figures/Fig4_trends/run_fig4.sh", "Fig 4 and its Mann-Whitney / Kruskal-Wallis statistics", timeout=1800),
    S("s4_fig5", "4", "bash figures/Fig5_proteome_scaling/run_fig5.sh", "Fig 5 power-law fits, level differences, clade-bootstrap CIs, lineage statistics", timeout=7200),
    S("s5_clusters", "5", "bash 5_umap/run_clusters.sh", "Dark / metazoan-protein / vertebrate-DNA cluster selection and Fig 7", timeout=3600),
    # ---------------------------------------------------------------- Stage 6: dark-cluster tables and Fig 8 from the archived CLEAN and Horizyn results
    S("s6_tables", "6", "bash 6_dark_cluster_annotation/run_tables.sh", "Horizyn family table, CLEAN summary, Tables S1 and S2, Fig 7 with the Fig 8 clusters", timeout=3600),
    S("s6_fig8", "6", "bash figures/Fig8_reaction_profiles/run_fig8.sh", "Fig 8: Horizyn reaction profiles of the rule-chosen dark clusters", timeout=1800),
    # ---------------------------------------------------------------- Stage 7: motifs from the archived STREME / FIMO results
    S("s7_no_motif_umap", "7", "bash 7_motif_analysis/run_no_motif_umap.sh", "Fig 11: enzymes with no discriminative motif on the map", timeout=1800),
    S("s7_motif_pick", "7", "bash 7_motif_analysis/run_motif_structure_pick.sh", "Motif structure pick by rule, Table 2", timeout=3600),
    # ---------------------------------------------------------------- the numbers and figure index
    S("s9_key_numbers", "9", "python workflow/figure_index.py && python workflow/key_numbers.py", "key_numbers.md (every number quoted in the text) and the figure index", timeout=900),

    # ---------------------------------------------------------------- tier B: GPU
    S("b_deploy", "2", "CLUSTER_WHAT=sam bash 1_clustering/run_mmseqs_clusters.sh && bash 2_training/run_deploy_ensemble.sh", "Deployment ensemble re-trained from the labels: thresholds, out-of-fold MCC, out-of-fold scores", tier="B", check="deployment", timeout=7200),
]

BY_ID = {s["id"]: s for s in STEPS}
