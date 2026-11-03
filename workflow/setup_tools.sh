#!/bin/bash
# Check (and, on request, install) the tools the pipeline needs that `make env` does not provide. Details and the reasons: EXTERNAL_TOOLS.md.
#   bash workflow/setup_tools.sh                  check everything, print what is missing and how to fix it (changes nothing)
#   bash workflow/setup_tools.sh --install foldseek|clean|horizyn     install one of them (network: clones / downloads, see EXTERNAL_TOOLS.md)
# ColabFold and ChimeraX are not installed by this script (large vendor installers): their steps are in EXTERNAL_TOOLS.md.
# Where the pipeline looks for each tool: refresh_config.json ("tools", "annotation"); an environment variable of the same name in capitals overrides it.
cd "$(dirname "$0")/.." || exit 1
REPO="$PWD"
cfg() { python3 -c "import json,sys;c=json.load(open('refresh_config.json'));print(c.get('annotation',{}).get(sys.argv[1]) or c.get('tools',{}).get(sys.argv[1]) or '')" "$1"; }
ok=0; miss=0
say() { if [ "$1" = ok ]; then echo "  [ok]      $2"; ok=$((ok + 1)); else echo "  [MISSING] $2"; echo "            -> $3"; miss=$((miss + 1)); fi; }
have_env() { conda env list 2> /dev/null | awk '{print $1}' | grep -qx "$1"; }

install_foldseek() { make tools; }
install_clean() {
    local dir="${CLEAN_DIR:-$REPO/tools/CLEAN}"
    [ -d "$dir/.git" ] || git clone https://github.com/tttianhao/CLEAN.git "$dir"
    have_env clean || conda env create -f 6_dark_cluster_annotation/clean/environment_clean.yml
    echo "CLEAN cloned to $dir. Set annotation.clean_app_dir in refresh_config.json to $dir/app. The pretrained CLEAN weights (app/data/pretrained) come with the repository;"
    echo "ESM-1b (about 2.5 GB) is fetched by CLEAN itself on its first run. See EXTERNAL_TOOLS.md for the files this pipeline expects in app/."
}
install_horizyn() {
    local dir="${HORIZYN_DIR:-$REPO/tools/horizyn}"
    [ -d "$dir/.git" ] || git clone https://github.com/dayhofflabs/horizyn.git "$dir"
    have_env horizyn || conda env create -f 6_dark_cluster_annotation/horizyn/environment_horizyn.yml
    mkdir -p "$dir/checkpoints"
    for f in horizyn_v1_0_dev.ckpt horizyn_v1_0_inf.ckpt; do        # the two official checkpoints (about 201 MB each), Zenodo record 20348783
        [ -s "$dir/checkpoints/$f" ] || curl -L -C - --retry 8 --retry-delay 10 -o "$dir/checkpoints/$f" "https://zenodo.org/api/records/20348783/files/$f/content"
    done
    (cd "$dir/checkpoints" && md5sum *.ckpt); echo "expected: dev=5b1f938f8b0a82fbe91892a3b4e2bf2c  inf=cf6775b775287462099ae0681485a6bc"
    echo "Set annotation.horizyn_dir in refresh_config.json to $dir (the folder with checkpoints/, configs/, horizyn/)."
}
if [ "$1" = "--install" ]; then
    case "$2" in foldseek) install_foldseek ;; clean) install_clean ;; horizyn) install_horizyn ;; *) echo "install: foldseek | clean | horizyn" >&2; exit 2 ;; esac
    exit $?
fi

echo "Tools of the pipeline (refresh_config.json, EXTERNAL_TOOLS.md)"
for t in fimo streme mmseqs hmmbuild mafft; do command -v $t > /dev/null && say ok "$t (bioconda, from environment.yml)" || say miss "$t" "make env   (and conda activate smstpred)"; done
fs="${FOLDSEEK_BIN:-$(cfg foldseek_bin)}"; [ -x "$fs" ] || fs=""; [ -z "$fs" ] && fs="$(command -v foldseek || true)"; [ -z "$fs" ] && [ -x "$REPO/tools/foldseek/bin/foldseek" ] && fs="$REPO/tools/foldseek/bin/foldseek"
[ -n "$fs" ] && say ok "foldseek: $fs" || say miss "foldseek (bake-off baseline, SaProt embeddings)" "bash workflow/setup_tools.sh --install foldseek      (= make tools)"
cb="${COLABFOLD_BATCH:-$(cfg colabfold_batch)}"
[ -n "$cb" ] && [ -x "$cb" ] && say ok "colabfold_batch: $cb" || say miss "colabfold_batch (structures AlphaFold DB cannot supply)" "install localcolabfold (EXTERNAL_TOOLS.md, section ColabFold), then set tools.colabfold_batch in refresh_config.json"
cx="${CHIMERAX:-$(cfg chimerax)}"
[ -n "$cx" ] && [ -x "$cx" ] && say ok "ChimeraX: $cx" || say miss "ChimeraX (Figs 9 and 10 panels)" "install ChimeraX (EXTERNAL_TOOLS.md, section ChimeraX), then set tools.chimerax in refresh_config.json"
ca="${CLEAN_APP_DIR:-$(cfg clean_app_dir)}"
{ [ -n "$ca" ] && [ -f "$ca/CLEAN_infer_fasta.py" ]; } && say ok "CLEAN checkout: $ca" || say miss "CLEAN checkout (Stage 6)" "bash workflow/setup_tools.sh --install clean, then set annotation.clean_app_dir"
have_env "$(cfg clean_env || echo clean)" && say ok "conda environment for CLEAN" || say miss "conda environment 'clean'" "conda env create -f 6_dark_cluster_annotation/clean/environment_clean.yml"
hz="${HORIZYN_DIR:-$(cfg horizyn_dir)}"
{ [ -n "$hz" ] && [ -s "$hz/checkpoints/horizyn_v1_0_inf.ckpt" ]; } && say ok "Horizyn checkout with checkpoints: $hz" || say miss "Horizyn checkout / checkpoints (Stage 6)" "bash workflow/setup_tools.sh --install horizyn, then set annotation.horizyn_dir"
have_env "$(cfg horizyn_env || echo horizyn)" && say ok "conda environment for Horizyn" || say miss "conda environment 'horizyn'" "conda env create -f 6_dark_cluster_annotation/horizyn/environment_horizyn.yml"
echo; echo "$ok ok, $miss missing"
[ "$miss" -eq 0 ]
