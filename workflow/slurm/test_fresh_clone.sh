#!/bin/bash
#SBATCH --job-name=smstpred_fresh_clone
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=fresh_clone_%j.out

# Fresh-clone test: export the repository as a commit of the current work would contain it (tracked + new files, nothing ignored), build the
# conda environment from environment.yml, put the Zenodo data in place, and run `make verify` (plus `make verify-gpu` with GPU=1).
# Submit from the repository root:   sbatch workflow/slurm/test_fresh_clone.sh
#   DATA_SRC   zenodo:10.5281/zenodo.23655998 (default), or a folder holding the record's files (stands in for an unpublished record)
#   TEST_DIR   where the clone is built (default: ../_fresh_clone_test next to the repository)
#   ENV_NAME   conda environment to create (default: smstpred_fresh; an existing one of that name is replaced)
#   GPU=1      also run `make verify-gpu`
# Results (verify/report.json, logs) are copied back into the repository's verify/ folder.
set -eo pipefail          # no -u: some conda activation hooks read unset variables
SRC="${SLURM_SUBMIT_DIR:-$PWD}"
DATA_SRC="${DATA_SRC:-zenodo:10.5281/zenodo.23655998}"
T="${TEST_DIR:-$SRC/../_fresh_clone_test}"
ENVN="${ENV_NAME:-smstpred_fresh}"
[ -f "$SRC/environment.yml" ] || { echo "submit from the repository root (no environment.yml in $SRC)" >&2; exit 1; }

rm -rf "$T"; mkdir -p "$T"
(cd "$SRC" && git ls-files --cached --others --exclude-standard -z) | rsync -a --from0 --files-from=- --ignore-missing-args "$SRC/" "$T/"
cd "$T"

eval "$(conda shell.bash hook)"
conda env remove -y -n "$ENVN" >/dev/null 2>&1 || true
conda env create -f environment.yml -n "$ENVN"
conda activate "$ENVN"
python --version; python -c "import numpy,pandas,scipy,sklearn,torch,matplotlib;print('numpy',numpy.__version__,'pandas',pandas.__version__,'scipy',scipy.__version__,'sklearn',sklearn.__version__,'torch',torch.__version__,'mpl',matplotlib.__version__)"

case "$DATA_SRC" in zenodo:*) make data DATA_SRC="$DATA_SRC" ;; *) make data DATA_SRC="$DATA_SRC" LINK=--link ;; esac
make verify || true
if [ "${GPU:-0}" = 1 ]; then make verify-gpu || true; fi

mkdir -p "$SRC/verify"
rsync -a "$T/verify/" "$SRC/verify/"
echo "report copied to $SRC/verify/report.json"
