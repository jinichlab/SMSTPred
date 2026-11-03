#!/bin/bash
#SBATCH --job-name=dark_tables
#SBATCH --ntasks=1
#SBATCH --output=run_tables_output.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G

# Tables S1 and S2 for a release's dark clusters, from the CLEAN and Horizyn results of <release>/dark_annotation (see the README of this folder):
#   clean/dark53_maxsep.csv + dark53_pvalue.csv   (clean/run_clean_dark.sh)     horizyn_run/horizyn_full_topk.csv   (horizyn/run_horizyn_predict_full.sh)
#   RELEASE_DIR=<release> [DARK_NAME=dark53] [EC_NAMES=<tab-separated EC, name file>] sbatch 6_dark_cluster_annotation/run_tables.sh
# Then Fig 8 (the clusters where CLEAN and Horizyn agree: figures/Fig8_reaction_profiles/select_panels.py) and Fig 7 again with those clusters outlined.
# Writes <release>/fig8/, <release>/dark_annotation/{clean/clean_summary.csv, horizyn/protein_clusters_family_final.csv, horizyn/avg_top1_score.csv, table_dark_clusters.csv, tables/table_s1.*, table_s2.*}.
here="$(dirname "$0")"
[ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
cd "$here" || exit 1
HERE="$PWD"
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
set -e
RELEASE_DIR="${RELEASE_DIR:-$(python ../workflow/release_dir.py)}"      # unset: the newest release
A="$RELEASE_DIR/dark_annotation"
NAME="${DARK_NAME:-$(cat "$A/dark_name.txt")}"                       # dark_name.txt: written by make_dark_inputs.py
python horizyn/final_family_table.py --topk "$A/horizyn_run/horizyn_full_topk.csv" --out "$A/horizyn"
# enzyme names: EC_NAMES (a tab-separated EC<TAB>name file, or an enzyme.dat) if set, else the ExPASy ENZYME file that 0_uniprot_pull/fetch_reference_tables.py put in the release
[ -z "${EC_NAMES:-}" ] && [ -f "$RELEASE_DIR/reference_tables/enzyme.dat" ] && EC_NAMES="$RELEASE_DIR/reference_tables/enzyme.dat"
python clean/summarize_clean_clusters.py --maxsep "$A/clean/${NAME}_maxsep.csv" --pvalue "$A/clean/${NAME}_pvalue.csv" --members "$A/${NAME}_members.csv" --out "$A/clean" ${EC_NAMES:+--ec-names "$EC_NAMES"}
python context_notes.py --release "$RELEASE_DIR" > "$A/context_notes.log" 2>&1 || echo "context notes skipped (UniProt not reachable?): the table is built without pathway notes"      # the automatic "-> pathway" notes (UniProt; network)
python build_table2.py --release-dir "$A"
python ../figures/TableS1_actinomycetota_dark_clusters/build_table.py --src "$A/table_dark_clusters.csv" --out "$A/tables"
python ../figures/TableS2_plant_dark_clusters/build_table.py --src "$A/table_dark_clusters.csv" --out "$A/tables"
python ../figures/Fig8_reaction_profiles/make_fig8_reaction_profiles.py --release "$RELEASE_DIR"
python ../figures/Fig7_clusters/build_cluster_grids.py --release "$RELEASE_DIR"       # again, now outlining the clusters of Fig 8
python ../figures/Fig7_clusters/make_fig7_clusters.py --dir "$RELEASE_DIR/umap/fig7"
