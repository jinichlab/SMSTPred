#!/bin/bash
# Refresh the whole analysis for the newest UniProt release as ONE SLURM job chain, with no input in between:
#
#   pull -> bin -> cluster (+ the 16 SAM train / validation / test splits) and cluster_up -> up_splits -> structures -> embeddings
#        -> bakeoff (16-task GPU array) and deploy (the deployment ensemble) -> report ;  sam_vs_up (16-task GPU array) -> sam_vs_up_report
#   then the analyses of the release, on the deployed ensemble's calls:
#        proteome embeddings -> predict -> tree (Fig 3), fig4, fig5, umap -> clusters (Fig 7) -> dark_inputs -> clean + horizyn -> annotate (Tables S1/S2, Fig 8, Fig 7 again)
#        predict -> domain_table (Table S3) -> dedup -> streme (4 classes side by side) -> fimo (motif specificity matrix) -> motif_structures (representing structure of every motif)
#
#   bash workflow/run_refresh.sh                    # everything, from a fresh UniProt pull
#   bash workflow/run_refresh.sh --from bin         # start at a stage (the earlier ones are assumed done for the newest release)
#   bash workflow/run_refresh.sh --from bakeoff --to report
#   bash workflow/run_refresh.sh --dry-run          # print the jobs and their dependencies, submit nothing
#   make refresh   /   make refresh-dry-run         (REFRESH_ARGS="--from bin" to pass options)
#
# Every stage is its own job (the scripts named below) and starts only when the jobs it needs have FINISHED WITHOUT ERROR (--dependency=afterok): a failed stage stops the chain and
# everything after it is cancelled (--kill-on-invalid-dep=yes), so nothing downstream runs on incomplete input.
#
#   stage            what runs                                                              waits for
#   pull             0_uniprot_pull/run_pull.sh                      (network)               -
#   bin              0_uniprot_pull/run_bin.sh                                               pull
#   cluster          1_clustering/run_mmseqs_clusters.sh  (SAM; + tvt_splits.py)             bin
#   cluster_up       1_clustering/run_mmseqs_clusters.sh  (CLUSTER_WHAT=up, 32 cpus, 100 G)  bin
#   up_splits        1_clustering/up_splits.py                                               cluster, cluster_up
#   structures       embeddings/fetch_structures.py                  (network)               bin
#   colabfold        embeddings/predict_structures.sh   (ColabFold for the proteins without an AlphaFold DB model; skipped in seconds when there are none; ingests into the structure store)  structures
#   embeddings       embeddings/run_embeddings.sh x3 (GPU): sam (ProtT5, ProstT5), structures (SaProt, ESM-IF), up-benchmark (ProtT5 of the UP set)    bin / structures / bin
#   bakeoff          2_training/run_bakeoff.sh          (GPU array, 16 splits)               cluster, embeddings (sam, structures)
#   deploy           2_training/run_deploy_ensemble.sh  (6-member ProtT5 ensemble)             cluster, embeddings (sam)
#   report           2_training/run_bakeoff_report.sh   (Fig 2A, Table 1, tests)             bakeoff
#   sam_vs_up        2_training/run_sam_vs_up.sh        (GPU array, 16 splits)               bakeoff, up_splits, embeddings (up-benchmark)
#   sam_vs_up_report 2_training/run_sam_vs_up_report.sh (paired tests of MTase-trained vs UP-trained)   sam_vs_up
#   publish_classifier workflow/publish_classifier.py (the newest ensemble -> ready_classifier/model/ + LATEST.json: the ready classifier always uses the newest trained model; the chain never pushes)   deploy
#   fig2             2_training/run_fig2.sh             (Fig 2: A representations and tools + B MTase-trained vs UP-trained, one figure)   report, sam_vs_up_report
#   ref_tables       0_uniprot_pull/run_reference_tables.sh (ExPASy ENZYME + Rhea EC links; network)    pull
#   embed_proteomes  embeddings/run_embeddings.sh EMBED_WHAT=proteomes (ProtT5 of the reference proteomes' MTase-domain proteins)   pull
#   predict          3_predictions/run_predict_release.sh (+ combine_predictions.py)         deploy, embed_proteomes
#   tree             4_phylo_trends/run_tree.sh         (Fig 3, no iTOL)                     predict
#   fig4, fig5       figures/Fig4_trends/run_fig4.sh, figures/Fig5_proteome_scaling/run_fig5.sh   predict
#   umap             5_umap/run_umap.sh                 (global UMAP + HDBSCAN; an existing map of the same proteins is REUSED: UMAP is not reproducible, UMAP_RECOMPUTE=1 forces a new one)   predict, embed_proteomes
#   clusters         5_umap/run_clusters.sh             (dark / metazoan-protein / vertebrate-DNA clusters, Fig 7)   umap
#   dark_inputs      6_dark_cluster_annotation/run_dark_inputs.sh   (dark proteins -> FASTA + members)   clusters
#   clean, horizyn   6_dark_cluster_annotation/run_clean_release.sh, run_horizyn_release.sh  (GPU)       dark_inputs
#   annotate         6_dark_cluster_annotation/run_tables.sh (CLEAN summary, Horizyn families, Tables S1/S2, Fig 8, Fig 7 with Fig 8's clusters outlined)   clean, horizyn, ref_tables
#   domain_table     7_motif_analysis/run_domain_table.sh (corpus, InterPro enrichment, Table S3)      predict
#   dedup, streme    7_motif_analysis/run_dedup.sh; run_streme.sh (one prep job, then one per class)    domain_table / dedup
#   fimo             7_motif_analysis/run_fimo_cross.sh (32 motifs x every protein -> specificity matrix, the 2-fold cut)   streme
#   motif_status     7_motif_analysis/run_fimo_all.sh -> run_no_motif_umap.sh (each protein against its own class's motifs; Fig 11, the no-motif UMAP)   streme, umap
#   motif_structures 7_motif_analysis/run_motif_structure_search.sh -> run_motif_structure_pick.sh (the structure of every motif >= 2-fold, by rule, with SAM + substrate; Table 2; network)   fimo
#   figure_index     workflow/figure_index.py + key_numbers.py (FIGURES.md / FIGURES.html, article_figures/, key_numbers.md of the release: what to review and update in the draft)   everything above
#   fig9_10          figures/Fig9_10_structures/run_fig9_10.sh (four closest / four furthest motifs by rule, rendered with ChimeraX)   motif_structures
#
# NOT automatic, and said so rather than hidden:
#   * proteins that AlphaFold DB cannot supply a structure for need ColabFold (not installed by this repository: embeddings/README.md, "Installing ColabFold"). The structures job writes
#     <release>/structures/to_predict.fasta; if it is not empty, run `sbatch embeddings/predict_structures.sh` and `python embeddings/fetch_structures.py --ingest-predictions ...`
#     before `--from embeddings`. The bake-off prints a WARNING with the number of benchmark proteins that lack a vector and are left out.
#   * Foldseek (`make tools`), the conda environments of `make env` and `make embeddings-env`.
# Run it from an ACTIVATED environment (conda activate smstpred): the jobs inherit it. REFERENCE_STRUCTURES=<folder of archived <Entry>.pdb> lets the structures stage reuse them where
# the sequence still matches. UNIPROT_PULL_DIR / SMSTPRED_DATA / ARCHIVED_EMBEDDINGS ... are the variables of the stage scripts (their headers).
# Settings of the splits, the MLPs, the threshold and the deployed model (control): refresh_config.json (copied into <release>/benchmark/config_used.json by the bake-off).
# Not in the chain: Figs 1 and S1 (dropped from the paper). The protein-cellularity violin is Fig 4A (figures/Fig4_trends/build_fig4_panels.py); 4_phylo_trends/protein_cellularity_violin.py is its archived original.
# lineage_statistics.py runs at the end of run_fig5.sh; the no-motif UMAP (Fig 11) is stage motif_status. The CLEAN and Horizyn checkouts are outside this repository: refresh_config.json "annotation".
# The Fig 8 clusters are chosen by a rule (figures/Fig8_reaction_profiles/select_panels.py), the panel titles come from panel_titles.csv: a reaction not listed there gets "<CLEAN enzyme name> (RHEA:id)".
set -e
cd "$(dirname "$0")/.."

STAGES=(pull bin cluster cluster_up up_splits structures colabfold embeddings bakeoff deploy report sam_vs_up sam_vs_up_report fig2 publish_classifier
        ref_tables embed_proteomes predict tree fig4 fig5 umap clusters dark_inputs clean horizyn annotate domain_table dedup streme fimo motif_status motif_structures fig9_10 figure_index)
FROM=pull; TO=figure_index; DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        --to) TO="$2"; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) sed -n '2,/^set -e/p' "$0" | sed 's/^# \{0,1\}//' | sed '$d'; exit 0 ;;
        *) echo "unknown option $1 (--from STAGE --to STAGE --dry-run)" >&2; exit 2 ;;
    esac
done
idx() { for i in "${!STAGES[@]}"; do [ "${STAGES[$i]}" = "$1" ] && { echo "$i"; return; }; done; echo "no such stage: $1 (${STAGES[*]})" >&2; exit 2; }
I_FROM=$(idx "$FROM"); I_TO=$(idx "$TO")
[ "$I_FROM" -le "$I_TO" ] || { echo "--from must not come after --to" >&2; exit 2; }
command -v sbatch > /dev/null || [ "$DRY" = 1 ] || { echo "sbatch not found: this chain needs SLURM (use --dry-run to see the plan)" >&2; exit 1; }

declare -A JOB
# submit NAME "DEP STAGES" cmd... : submit one job after the jobs of the stages in DEP that are part of this chain
submit() {
    local name="$1" deps="$2"; shift 2
    local dep="" d
    for d in $deps; do [ -n "${JOB[$d]:-}" ] && dep="$dep:${JOB[$d]}"; done
    local opt=(--parsable --export=ALL); [ -n "$dep" ] && opt+=(--dependency="afterok$dep" --kill-on-invalid-dep=yes)
    if [ "$DRY" = 1 ]; then
        JOB[$name]="<$name>"
        printf '%-16s %s%ssbatch %s %s\n' "$name" "${EMBED_WHAT:+EMBED_WHAT=$EMBED_WHAT }" "${CLUSTER_WHAT:+CLUSTER_WHAT=$CLUSTER_WHAT }" "${opt[*]}" "$*"
    else
        JOB[$name]=$(sbatch "${opt[@]}" "$@")
        printf '%-16s job %s  (after:%s)\n' "$name" "${JOB[$name]}" "${dep:- nothing}"
    fi
}

log() { echo "--output=$PWD/$1"; }
wanted() { local i; i=$(idx "$1"); [ "$i" -ge "$I_FROM" ] && [ "$i" -le "$I_TO" ]; }
echo "refresh chain: $FROM -> $TO$( [ "$DRY" = 1 ] && echo '   (dry run, nothing is submitted)')"
REFS="${REFERENCE_STRUCTURES:+--reuse $REFERENCE_STRUCTURES}"
wanted pull       && submit pull       ""                    0_uniprot_pull/run_pull.sh
wanted bin        && submit bin        "pull"                0_uniprot_pull/run_bin.sh
wanted cluster    && CLUSTER_WHAT=sam submit cluster "bin"   1_clustering/run_mmseqs_clusters.sh
wanted cluster_up && CLUSTER_WHAT=up  submit cluster_up "bin" --cpus-per-task=32 --mem=100G 1_clustering/run_mmseqs_clusters.sh
wanted up_splits  && submit up_splits  "cluster cluster_up"  --cpus-per-task=2 --mem=16G --job-name=up_splits --output=up_splits_output.log --chdir=1_clustering --wrap "python -u up_splits.py"
wanted structures && submit structures "bin"                 --cpus-per-task=4 --mem=16G --job-name=fetch_structures --output=fetch_structures_output.log --chdir=embeddings --wrap "for i in 1 2 3 4 5 6; do python -u fetch_structures.py $REFS && exit 0; echo \"retrying the structure fetch (\$i)\"; sleep 180; done; exit 1"
wanted colabfold  && submit colabfold  "structures"          $(log embeddings/predict_structures_output.log) embeddings/predict_structures.sh      # ColabFold for the proteins AlphaFold DB cannot supply (exits at once when none)
if wanted embeddings; then
    EMBED_WHAT=sam          submit embed_sam        "bin"        embeddings/run_embeddings.sh
    EMBED_WHAT=structures   submit embed_structures "colabfold" embeddings/run_embeddings.sh
    EMBED_WHAT=up-benchmark submit embed_up         "bin"        embeddings/run_embeddings.sh
fi
wanted bakeoff    && submit bakeoff    "cluster embed_sam embed_structures"   2_training/run_bakeoff.sh
wanted deploy     && submit deploy     "cluster embed_sam"                    2_training/run_deploy_ensemble.sh
wanted report     && submit report     "bakeoff"                    2_training/run_bakeoff_report.sh
wanted sam_vs_up  && submit sam_vs_up  "bakeoff up_splits embed_up"           2_training/run_sam_vs_up.sh
wanted sam_vs_up_report && submit sam_vs_up_report "sam_vs_up"                2_training/run_sam_vs_up_report.sh
wanted publish_classifier && submit publish_classifier "deploy" --cpus-per-task=1 --mem=2G $(log workflow/publish_classifier_output.log) --wrap "python workflow/publish_classifier.py"
wanted fig2             && submit fig2             "report sam_vs_up_report"   2_training/run_fig2.sh

# ---- the analyses of the release: proteome predictions -> Fig 3/4/5, UMAP + clusters (Fig 7), dark-cluster annotation (Tables S1/S2, Fig 8), motif analysis (Table S3, STREME, FIMO)
# These scripts find the release themselves when they run (the newest in $UNIPROT_PULL_DIR, or $RELEASE_DIR), so the chain can be submitted before the pull has made its folder.
wanted ref_tables       && submit ref_tables       "pull"                      $(log 0_uniprot_pull/reference_tables_output.log)     0_uniprot_pull/run_reference_tables.sh
wanted embed_proteomes  && EMBED_WHAT=proteomes submit embed_proteomes "pull"  embeddings/run_embeddings.sh
wanted predict          && submit predict          "deploy embed_proteomes"    $(log 3_predictions/predict_release_output.log)       3_predictions/run_predict_release.sh
wanted tree             && submit tree             "predict"                   $(log 4_phylo_trends/run_tree_output.log)             4_phylo_trends/run_tree.sh
wanted fig4             && submit fig4             "predict"                   $(log figures/Fig4_trends/run_fig4_output.log)        figures/Fig4_trends/run_fig4.sh
wanted fig5             && submit fig5             "predict"                   $(log figures/Fig5_proteome_scaling/run_fig5_output.log) figures/Fig5_proteome_scaling/run_fig5.sh
wanted umap             && submit umap             "predict embed_proteomes"   $(log 5_umap/run_umap_output.log)                     5_umap/run_umap.sh
wanted clusters         && submit clusters         "umap"                      $(log 5_umap/run_clusters_output.log)                 5_umap/run_clusters.sh
wanted dark_inputs      && submit dark_inputs      "clusters"                  $(log 6_dark_cluster_annotation/run_dark_inputs_output.log) 6_dark_cluster_annotation/run_dark_inputs.sh
wanted clean            && submit clean            "dark_inputs"               $(log 6_dark_cluster_annotation/clean_dark.log)       6_dark_cluster_annotation/run_clean_release.sh
wanted horizyn          && submit horizyn          "dark_inputs"               $(log 6_dark_cluster_annotation/horizyn_dark.log)     6_dark_cluster_annotation/run_horizyn_release.sh
wanted annotate         && submit annotate         "clean horizyn ref_tables"  $(log 6_dark_cluster_annotation/run_tables_output.log) 6_dark_cluster_annotation/run_tables.sh
wanted domain_table     && submit domain_table     "predict"                   $(log 7_motif_analysis/domain_table_output.log)       7_motif_analysis/run_domain_table.sh
wanted dedup            && MOTIF_DIR=auto submit dedup "domain_table"          $(log 7_motif_analysis/dedup_output.log)              7_motif_analysis/run_dedup.sh
if wanted streme; then   # STREME is single-threaded: one job prepares the inputs, then the four classes run side by side
    MOTIF_DIR=auto STREME_PREP=1 STREME_CLASSES="" submit streme_prep "dedup" $(log 7_motif_analysis/streme_prep_output.log) 7_motif_analysis/run_streme.sh
    for cls in Metabolite Protein RNA DNA; do
        MOTIF_DIR=auto STREME_PREP=0 STREME_CLASSES="$cls" submit "streme_$cls" "streme_prep" $(log 7_motif_analysis/streme_${cls}_output.log) 7_motif_analysis/run_streme.sh
    done
fi
wanted fimo             && MOTIF_DIR=auto submit fimo "streme_Metabolite streme_Protein streme_RNA streme_DNA" --cpus-per-task=4 --mem=16G $(log 7_motif_analysis/fimo_cross_output.log) 7_motif_analysis/run_fimo_cross.sh
if wanted motif_status; then    # every protein against its own class's motifs (motif_status_all.csv) -> the no-motif UMAP (Fig 11)
    MOTIF_DIR=auto submit motif_status "streme_Metabolite streme_Protein streme_RNA streme_DNA" --cpus-per-task=4 --mem=16G $(log 7_motif_analysis/fimo_all_output.log) 7_motif_analysis/run_fimo_all.sh
    MOTIF_DIR=auto submit no_motif_umap "motif_status umap" $(log 7_motif_analysis/no_motif_umap_output.log) 7_motif_analysis/run_no_motif_umap.sh
fi
if wanted motif_structures; then   # which carriers of each motif have a structure with SAM + substrate (RCSB/UniProt metadata), then the representing structure of every motif by rule (mmCIF download)
    MOTIF_DIR=auto submit motif_search "fimo"        --cpus-per-task=2 --mem=16G $(log 7_motif_analysis/motif_structure_search_output.log) 7_motif_analysis/run_motif_structure_search.sh
    MOTIF_DIR=auto submit motif_pick   "motif_search" --cpus-per-task=2 --mem=16G $(log 7_motif_analysis/motif_structure_pick_output.log)   7_motif_analysis/run_motif_structure_pick.sh
fi
wanted fig9_10          && submit fig9_10          "motif_pick"                $(log figures/Fig9_10_structures/run_fig9_10_output.log) figures/Fig9_10_structures/run_fig9_10.sh
wanted figure_index     && submit figure_index     "report sam_vs_up_report fig2 tree fig4 fig5 annotate no_motif_umap motif_pick fig9_10 domain_table" --cpus-per-task=1 --mem=2G $(log workflow/figure_index_output.log) --wrap "python workflow/figure_index.py && python workflow/key_numbers.py"
echo
[ "$DRY" = 1 ] || echo "follow with: squeue -u \$USER ; the results appear in <release>/benchmark/ (Fig 2A, Table 1, Fig S4, paired tests, loss curves) and <release>/sam_vs_up/"
