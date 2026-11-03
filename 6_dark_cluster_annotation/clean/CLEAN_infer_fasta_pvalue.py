"""CLEAN p-value inference for a named FASTA (data/inputs/<name>.fasta of the CLEAN checkout); the repo's copy of the checkout's CLEAN_infer_fasta_pvalue.py, which is hard-wired to dark33.
Run from the app/ directory of the CLEAN checkout (run_clean_dark.sh does): python CLEAN_infer_fasta_pvalue.py -d dark53. Writes results/inputs/<name>_pvalue.csv."""
import argparse

from CLEAN.infer import infer_pvalue
from CLEAN.utils import prepare_infer_fasta

ap = argparse.ArgumentParser()
ap.add_argument("-d", "--fasta_data", required=True)
name = "inputs/" + ap.parse_args().fasta_data
prepare_infer_fasta(name)            # dummy csv; the ESM embeddings are cached in data/esm_data by the maxsep run
infer_pvalue("split100", name, report_metrics=False, pretrained=True)
