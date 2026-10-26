#!/bin/bash
#SBATCH --job-name=fimo_cross
#SBATCH --ntasks=1
#SBATCH --output=fimo_cross_output.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00

# CROSS-SCAN: every class's 8 STREME motifs against EVERY corpus sequence.
#
# run_fimo_all.sh only ever scanned a class's motifs against that same class's
# members (a diagonal: streme_RNA x class_RNA.fasta). That answers "does a
# class-X protein carry class-X motifs?" (the 12.1%-motif-less result) but by
# construction it can NOT say how often a motif shows up OUTSIDE its class, and it
# never touched the 12,208 Mixed/Unassigned proteins. This run fills both gaps:
#
#   1. motif x class SPECIFICITY matrix -- the evidence behind calling these
#      motifs "discriminative" (so far that rests only on STREME discovery-time
#      enrichment, never verified by scanning). Rate-normalised, because the
#      classes differ ~6x in size (Other 41,941 vs DNA 6,822).
#   2. motif status for Mixed / Unassigned proteins (never scanned before).
#   3. a lookup to verify a structure's protein really carries its motif, for any
#      motif/protein pair regardless of class (see the structure-mapping audit).
#
# Same motif files and same p<1e-4 threshold as the diagonal, so coverage numbers
# stay directly comparable. --text skips the bulky cisml.xml/fimo.gff/html we
# never read (~38 MB/class); note q-values are not computed in text mode, so we
# threshold on p, exactly as the diagonal did.
#
# NOTE: corpus_all.fasta is NOT all proteins -- it is the 106,608 four-class
# members with '>Entry|Class' headers. We build the true all-protein FASTA
# (118,816, plain '>Entry' headers) from corpus.csv here.
set -e   # fail fast: don't cascade four broken FIMO runs off a bad prep

source ~/miniconda3/etc/profile.d/conda.sh
conda activate classifier_env

# SLURM copies this script to /var/spool/slurmd/..., so "$(dirname $0)" is NOT the
# script's real home -- use the submit dir (falling back to $0 for direct runs).
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
for need in corpus.csv streme_DNA/streme.txt streme_RNA/streme.txt \
            streme_Protein/streme.txt streme_Other/streme.txt; do
    [ -e "$need" ] || { echo "ERROR: $need not found in $PWD -- sbatch from motif_analysis/"; exit 1; }
done
echo "[cwd] $PWD"

# 1) build + sanitize the full target set (20 standard AAs; FIMO rejects U/O/B/Z)
python - <<'PY'
import re
import pandas as pd
c = pd.read_csv("corpus.csv", usecols=["Entry", "Sequence", "class"],
                keep_default_na=False, low_memory=False)
with open("corpus_every.clean.fasta", "w") as f:
    for e, s in zip(c["Entry"], c["Sequence"]):
        f.write(f">{e}\n{re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', 'X', str(s))}\n")
print(f"[prep] corpus_every.clean.fasta: {len(c):,} sequences "
      f"({(~c['class'].isin(['DNA','RNA','Protein','Other'])).sum():,} Mixed/Unassigned "
      f"the diagonal never scanned)")
PY

# 2) each class's 8 motifs vs ALL sequences
mkdir -p fimo_cross
for cls in DNA RNA Protein Other; do
    echo "=== FIMO cross ${cls}: 8 motifs vs $(grep -c '^>' corpus_every.clean.fasta) sequences ==="
    fimo --text --thresh 1e-4 --verbosity 1 \
        "streme_${cls}/streme.txt" corpus_every.clean.fasta > "fimo_cross/${cls}.tsv"
    echo "    $(( $(wc -l < "fimo_cross/${cls}.tsv") - 1 )) hits -> fimo_cross/${cls}.tsv"
done

# 3) specificity matrix + per-protein cross-class motif status
python fimo_cross_manifest.py

echo "fimo_cross done"
