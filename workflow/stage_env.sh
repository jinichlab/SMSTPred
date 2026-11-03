# Sourced by the stage scripts that run inside the refresh chain (workflow/run_refresh.sh) after the pull has made the release folder.
#   here="$(dirname "$0")"; [ "$(basename "$0")" = slurm_script ] && here="$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -n 1)")"
#   source "$here/../workflow/stage_env.sh"
# Gives the script: REPO (the repository), REL (the release folder: $RELEASE_DIR, else the newest in $UNIPROT_PULL_DIR), the active conda environment
# (the one already active, else ${SMSTPRED_ENV:-smstpred}), and refresh_config.json's "annotation" settings through `cfg KEY` (an environment variable of the same name in capitals wins).
REPO="$(cd "$here" && pwd)"; while [ ! -f "$REPO/workflow/release_dir.py" ] && [ "$REPO" != / ]; do REPO="$(dirname "$REPO")"; done      # the repository: the folder above the script that holds workflow/
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" = base ]; then
    if command -v conda > /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "${SMSTPRED_ENV:-smstpred}" || { echo "conda environment ${SMSTPRED_ENV:-smstpred} not found: run 'make env' first" >&2; exit 1; }
    fi
fi
REL="$(python "$REPO/workflow/release_dir.py")" || { echo "no release folder found (set RELEASE_DIR or UNIPROT_PULL_DIR)" >&2; exit 1; }
export REPO REL; export RELEASE_DIR="$REL"
echo "[stage] release $REL"
cfg() {   # cfg clean_app_dir -> $CLEAN_APP_DIR, else annotation.clean_app_dir (or tools.<key>) of refresh_config.json
    local up; up="$(echo "$1" | tr a-z A-Z)"
    if [ -n "${!up:-}" ]; then echo "${!up}"; else python -c "import json,sys;c=json.load(open(sys.argv[1]));print(c.get('annotation',{}).get(sys.argv[2]) or c.get('tools',{}).get(sys.argv[2]) or '')" "$REPO/refresh_config.json" "$1"; fi
}
