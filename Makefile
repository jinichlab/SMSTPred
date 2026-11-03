# One entry point for the whole analysis. `make help` lists the targets.
#
#   make env DATA_SRC=zenodo:10.5281/zenodo.23655998     # create the conda environment (once), then: conda activate smstpred
#   make data DATA_SRC=zenodo:10.5281/zenodo.23655998    # download + extract + verify the archived data (or DATA_SRC=/path LINK=--link)
#   make verify                                          # re-run every CPU step in a scratch copy and compare with the committed results

SHELL := /bin/bash
ENV_NAME ?= smstpred
PY ?= python
DATA_SRC ?=
LINK ?=
STEPS := $(PY) workflow/run_steps.py

.PHONY: help annotation-test pull pull-smoke pull-compare bin bin-test cluster splits up-splits structures bakeoff-test bakeoff bakeoff-report classical deploy sam-vs-up sam-vs-up-report tree refresh refresh-dry-run embeddings-env embeddings-check env data check tools steps stage1 stage2 stage3 stage4 stage5 stage6 stage7 figures tables verify verify-gpu report all clean-verify

help:
	@echo "env          create the conda environment ($(ENV_NAME)) from environment.yml"
	@echo "data         fetch the Zenodo data into data/   (DATA_SRC=zenodo:<doi> | DATA_SRC=/local/folder [LINK=--link])"
	@echo "check        verify the checksums of data/"
	@echo "pull         Stage 0 (network): pull the current UniProt SAM-MTase set, UP set and the 1,259 proteomes -> \$$UNIPROT_PULL_DIR or data/uniprot_pull/"
	@echo "pull-smoke   Stage 0 smoke test (3 proteomes, a few hundred rows)"
	@echo "pull-compare Stage 0: what changed between the newest pull and the archived data"
	@echo "bin          Stage 0: type the pulled proteins by substrate (Metabolite/Protein/RNA/DNA) -> labels_SAM.csv, labels_UP.csv, binning_report.txt"
	@echo "bin-test     Stage 0: unit tests of the substrate-typing rules (offline, seconds)"
	@echo "cluster      Stage 1: MMseqs2 30% clusters of the typed SAM sets + the 16 train/validation/test splits (after 'make bin'; the UP one: CLUSTER_WHAT=up)"
	@echo "splits       Stage 1: only the 16 cluster-disjoint splits (20% of the clusters test, 15% validation) of the newest release (seconds)"
	@echo "up-splits    Stage 1: the UP train/validation sets of the SAM-vs-UP comparison (UP clusters that contain a SAM test protein are removed; needs 'make cluster' with CLUSTER_WHAT=up)"
	@echo "bakeoff-test Stage 2: offline tests of the bake-off, the splits, the SAM-vs-UP and classical-ML code and the MCC threshold search (seconds)"
	@echo "bakeoff      Stage 2: submit the classifier bake-off, one SLURM array task per split (needs the embeddings, structures and splits of the release)"
	@echo "bakeoff-report  Stage 2: merge the 16 splits, statistics, Fig 2A, Table 1, paired tests, train/validation loss figures (and Fig S4 if 'make classical' has run)"
	@echo "classical    Stage 2: the classical-ML arm (XGBoost, RF, LR, 4 SVMs on the ProtT5 embeddings), one SLURM array task per split"
	@echo "deploy       Stage 2: the deployment model: a 6-member ProtT5 MLP ensemble from a 6-fold cluster partition of the deployment set, with honest out-of-fold predictions (SLURM, GPU, minutes)"
	@echo "sam-vs-up    Stage 2: SAM-trained vs UP-trained ProtT5 MLP on the same held-out SAM clusters, one SLURM array task per split (after 'make bakeoff')"
	@echo "sam-vs-up-report  Stage 2: the paired tests (Wilcoxon, Nadeau-Bengio) and figures of the SAM-vs-UP comparison"
	@echo "tree         Stage 4: Fig 3 for the newest release with no iTOL: combine proteome predictions, taxonomy tree + annotation files, draw the figure"
	@echo "refresh      the whole chain for the newest UniProt release as one SLURM job chain: pull, bin, cluster, structures, embeddings, bakeoff, report, then proteome predictions, Fig 3-5, UMAP + Fig 7, CLEAN + Horizyn + Tables S1/S2 + Fig 8, Table S3, STREME, FIMO (bash workflow/run_refresh.sh --help)"
	@echo "refresh-dry-run  print the job chain of 'make refresh' without submitting anything"
	@echo "structures   Embeddings: build the structure store (AlphaFold DB API, sequence-checked) and list what needs a ColabFold prediction"
	@echo "tools        download the Foldseek binary (only needed for the bake-off baseline and SaProt embeddings)"
	@echo "steps        list every step (id, tier, folder, command)"
	@echo "stage1 ... stage7   run that stage in place (CPU steps)"
	@echo "figures      run the figure and table builders in place"
	@echo "verify       run all CPU steps in a scratch copy and compare with the committed results -> verify/report.json"
	@echo "verify-gpu   re-train the deployment model, the MLP bake-off and the loss curves on a GPU and compare outcomes with stated tolerances"
	@echo "embeddings-env     create the two environments for re-computing embeddings (embeddings/)"
	@echo "embeddings-check   re-embed 20 proteins with each embedder and compare with the archived vectors (GPU)"
	@echo "report       render verify/report.json into the tables in REPRODUCING.md"
	@echo "all          data check, every stage, figures, then verify"

embeddings-env:   # the two environments of embeddings/ (only needed to re-compute embeddings)
	conda env create -f embeddings/environment_transformers.yml || conda env update -f embeddings/environment_transformers.yml
	conda env create -f embeddings/environment_esmif.yml || conda env update -f embeddings/environment_esmif.yml

embeddings-check:   # re-embed 20 proteins with every embedder and compare with the archived vectors (GPU)
	bash embeddings/run_spotcheck.sh

env:
	conda env create -f environment.yml -n $(ENV_NAME) || conda env update -f environment.yml -n $(ENV_NAME)
	@echo "now run: conda activate $(ENV_NAME)"

data:
	$(PY) workflow/fetch_data.py --src "$(DATA_SRC)" $(LINK)

check:
	$(PY) workflow/fetch_data.py --check

tools:
	mkdir -p tools && cd tools && wget -q https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && tar xzf foldseek-linux-avx2.tar.gz && rm foldseek-linux-avx2.tar.gz
	@echo "Foldseek is in tools/foldseek/bin (add it to PATH)"

pull:   # under SLURM use `sbatch 0_uniprot_pull/run_pull.sh` (same script)
	bash 0_uniprot_pull/run_pull.sh

pull-smoke:
	cd 0_uniprot_pull && $(PY) pull_uniprot.py all --limit 3

pull-compare:
	cd 0_uniprot_pull && $(PY) compare_to_previous.py

bin:    # under SLURM use `sbatch 0_uniprot_pull/run_bin.sh` (same script)
	bash 0_uniprot_pull/run_bin.sh

bin-test:
	cd 0_uniprot_pull && $(PY) test_substrate_binning.py

cluster:   # under SLURM use `sbatch 1_clustering/run_mmseqs_clusters.sh` (same script)
	bash 1_clustering/run_mmseqs_clusters.sh

splits:
	cd 1_clustering && $(PY) tvt_splits.py

up-splits:
	cd 1_clustering && $(PY) up_splits.py

bakeoff-test:
	cd 1_clustering && $(PY) test_tvt_splits.py && $(PY) test_up_splits.py
	cd 2_training && $(PY) test_bakeoff.py && $(PY) test_sam_vs_up.py && $(PY) test_deploy_ensemble.py && $(PY) test_threshold_search.py

bakeoff:   # a 16-task SLURM array (one split per task), about 10 minutes each; then `make bakeoff-report`
	sbatch 2_training/run_bakeoff.sh

bakeoff-report:   # under SLURM use `sbatch 2_training/run_bakeoff_report.sh` (same script)
	bash 2_training/run_bakeoff_report.sh

classical:   # a 16-task CPU array; the SVMs dominate the time
	sbatch 2_training/run_classical.sh

deploy:   # about 5 minutes on one GPU; writes <release>/deployment/ (DEPLOY_DIR=<release>/deployment for 3_predictions/predict.py)
	sbatch 2_training/run_deploy_ensemble.sh

sam-vs-up:   # a 16-task GPU array, about 5-10 minutes per split
	sbatch 2_training/run_sam_vs_up.sh

sam-vs-up-report:   # under SLURM use `sbatch 2_training/run_sam_vs_up_report.sh` (same script)
	bash 2_training/run_sam_vs_up_report.sh

tree:   # under SLURM use `sbatch 4_phylo_trends/run_tree.sh` (same script)
	bash 4_phylo_trends/run_tree.sh

annotation-test:   # the checks of the stages after the predictions: cluster selection, dark-protein inputs, CLEAN summary, Horizyn family table, Fig 8 rule
	cd 5_umap && $(PY) test_cluster_selection.py
	cd 6_dark_cluster_annotation && $(PY) test_make_dark_inputs.py
	cd 6_dark_cluster_annotation && $(PY) test_rhea_text.py
	cd 6_dark_cluster_annotation/clean && $(PY) test_summarize_clean.py
	cd 6_dark_cluster_annotation/horizyn && $(PY) test_final_family_table.py
	cd figures/Fig8_reaction_profiles && $(PY) test_select_panels.py
	cd 7_motif_analysis && $(PY) test_motif_structure_search.py && $(PY) test_motif_structure_pick.py

refresh:
	bash workflow/run_refresh.sh $(REFRESH_ARGS)

refresh-dry-run:
	bash workflow/run_refresh.sh --dry-run $(REFRESH_ARGS)

structures:   # REUSE=<folder of archived <Entry>.pdb> to reuse them where the sequence still matches
	cd embeddings && $(PY) fetch_structures.py $(if $(REUSE),--reuse $(REUSE),)

steps:
	$(STEPS) --list

stage1 stage2 stage3 stage4 stage5 stage6 stage7:
	$(STEPS) --mode run --stage $(subst stage,,$@)

figures tables:
	$(STEPS) --mode run --stage figures

verify:
	$(STEPS) --mode verify --tier A
	$(PY) workflow/make_report.py

verify-gpu:   # re-trains the deployment model and the MLP bake-off, re-predicts, re-runs the loss-curve sweeps; needs a GPU (uses sbatch if present)
	$(STEPS) --mode verify --tier B
	$(PY) workflow/make_report.py

report:
	$(PY) workflow/make_report.py

all: check stage1 stage2 stage3 stage4 stage5 stage6 stage7 figures verify

clean-verify:
	rm -rf .verify_work verify/logs
