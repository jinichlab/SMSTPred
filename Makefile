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

.PHONY: help embeddings-env embeddings-check env data check tools steps stage1 stage2 stage3 stage4 stage5 stage6 stage7 figures tables verify verify-gpu report all clean-verify

help:
	@echo "env          create the conda environment ($(ENV_NAME)) from environment.yml"
	@echo "data         fetch the Zenodo data into data/   (DATA_SRC=zenodo:<doi> | DATA_SRC=/local/folder [LINK=--link])"
	@echo "check        verify the checksums of data/"
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
