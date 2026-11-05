#!/bin/bash
# Check and install the tools the pipeline needs that `make env` does not provide. Details and the reasons: EXTERNAL_TOOLS.md.
#   bash workflow/setup_tools.sh                       check everything, print what is missing and how to fix it (changes nothing)
#   bash workflow/setup_tools.sh --install all         install every tool that can be installed without a licence step (foldseek, clean, horizyn, colabfold)
#   bash workflow/setup_tools.sh --install foldseek|clean|horizyn|colabfold
#   bash workflow/setup_tools.sh --install chimerax --chimerax-package <ChimeraX .deb | .tar.gz | https URL>
#       ChimeraX is the one tool a script cannot fetch by itself: UCSF hands out the package only after the licence form on https://www.cgl.ucsf.edu/chimerax/download.html
#       has been accepted. Accept it there (free for non-commercial use) and give this script the file or the link you were given; everything else is automatic.
#   bash workflow/setup_tools.sh --adopt KEY=PATH ...  record a tool that is already installed elsewhere (e.g. --adopt clean_app_dir=/opt/CLEAN/app)
# Everything is installed under <repo>/tools/ (override with TOOLS_DIR) and recorded in refresh_config.local.json (not committed), which the pipeline reads through
# workflow/cfgget.py (order: environment variable, refresh_config.local.json, refresh_config.json, <repo>/tools/). Needs outbound internet and conda.
cd "$(dirname "$0")/.." || exit 1
REPO="$PWD"
TOOLS="${TOOLS_DIR:-$REPO/tools}"
get() { python3 workflow/cfgget.py "$1"; }
record() {   # record KEY VALUE [section]  -> refresh_config.local.json
    python3 - "$1" "$2" "${3:-tools}" << 'PY'
import json, sys
from pathlib import Path
p = Path("refresh_config.local.json")
c = json.loads(p.read_text()) if p.exists() else {}
c.setdefault(sys.argv[3], {})[sys.argv[1]] = sys.argv[2]
p.write_text(json.dumps(c, indent=1) + "\n")
PY
}
ok=0; miss=0
say() { if [ "$1" = ok ]; then echo "  [ok]      $2"; ok=$((ok + 1)); else echo "  [MISSING] $2"; echo "            -> $3"; miss=$((miss + 1)); fi; }
have_env() { conda env list 2> /dev/null | awk '{print $1}' | grep -qx "$1"; }
need() { command -v "$1" > /dev/null || { echo "ERROR: $1 is required for this step" >&2; exit 1; }; }
fetch() { curl -fL --retry 8 --retry-delay 10 -C - -o "$2" "$1"; }

install_foldseek() {
    need curl; mkdir -p "$TOOLS/foldseek" && (cd "$TOOLS" && fetch https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz foldseek.tgz && tar xzf foldseek.tgz && rm foldseek.tgz) || exit 1
    "$TOOLS/foldseek/bin/foldseek" version && record foldseek_bin "$TOOLS/foldseek/bin/foldseek"
}
install_clean() {
    need git; need conda
    local dir="$TOOLS/CLEAN"
    [ -d "$dir/.git" ] || git clone https://github.com/tttianhao/CLEAN.git "$dir" || exit 1
    local ce; ce="$(get clean_env)"
    have_env "$ce" || { conda env create -n "$ce" -f 6_dark_cluster_annotation/clean/environment_clean.yml && conda run -n "$ce" pip install --no-deps -r 6_dark_cluster_annotation/clean/requirements_clean.txt; } || { conda env remove -y -n "$ce" > /dev/null 2>&1; exit 1; }
    [ -f "$dir/app/CLEAN_infer_fasta.py" ] || { echo "CLEAN checkout has no app/CLEAN_infer_fasta.py: the repository layout changed" >&2; exit 1; }
    # CLEAN's own install steps (its README, section 1.2/1.3): build the CLEAN package into the environment, vendor fair-esm (the pipeline puts app/esm on PYTHONPATH:
    # the environment's esm 3.x has another API), make the embedding folder, and fetch the pretrained weights from the link in the README (Google Drive, via gdown)
    (cd "$dir/app" && conda run -n "$ce" python build.py install) || exit 1
    [ -d "$dir/app/esm/.git" ] || { git clone https://github.com/facebookresearch/esm.git "$dir/app/esm" && git -C "$dir/app/esm" checkout -q 2b369911bb5b4b0dda914521b9475cad1656b2ac; } || exit 1
    mkdir -p "$dir/app/data/esm_data"
    if [ ! -s "$dir/app/data/pretrained/100.pt" ]; then
        mkdir -p "$dir/app/data/pretrained" "$dir/dl"
        conda run -n "$ce" gdown --fuzzy "https://drive.google.com/file/d/1kwYd4VtzYuMvJMWXy6Vks91DSUAOcKpZ/view" -O "$dir/dl/pretrained.zip" || { echo "could not fetch CLEAN's pretrained weights from Google Drive (quota?): download the file named in CLEAN's README by hand into $dir/dl/pretrained.zip and rerun" >&2; exit 1; }
        (cd "$dir/dl" && unzip -oq pretrained.zip) || exit 1
        f="$(find "$dir/dl" -name 100.pt | head -n 1)"; [ -n "$f" ] || { echo "100.pt not found in the CLEAN weights archive" >&2; exit 1; }
        cp "$(dirname "$f")"/*.pt "$(dirname "$f")"/*.pth "$(dirname "$f")"/*.pkl "$dir/app/data/pretrained/" 2> /dev/null
    fi
    [ -s "$dir/app/data/pretrained/100.pt" ] && [ -s "$dir/app/data/pretrained/70.pt" ] || { echo "CLEAN pretrained weights missing in $dir/app/data/pretrained" >&2; exit 1; }
    record clean_app_dir "$dir/app" annotation; record clean_env "$ce" annotation
    echo "CLEAN installed: $dir/app (conda environment '$ce')"
}
install_horizyn() {
    need git; need conda; need curl
    local dir="$TOOLS/horizyn"
    [ -d "$dir/.git" ] || git clone https://github.com/dayhofflabs/horizyn.git "$dir" || exit 1
    local he; he="$(get horizyn_env)"
    have_env "$he" || { conda env create -n "$he" -f 6_dark_cluster_annotation/horizyn/environment_horizyn.yml && conda run -n "$he" pip install --no-deps -r 6_dark_cluster_annotation/horizyn/requirements_horizyn.txt && conda run -n "$he" pip install --no-deps -e "$dir"; } || { conda env remove -y -n "$he" > /dev/null 2>&1; exit 1; }
    mkdir -p "$dir/checkpoints"
    for f in horizyn_v1_0_dev.ckpt horizyn_v1_0_inf.ckpt; do        # the two official checkpoints (about 201 MB each), Zenodo record 20348783
        [ -s "$dir/checkpoints/$f" ] || fetch "https://zenodo.org/api/records/20348783/files/$f/content" "$dir/checkpoints/$f" || exit 1
    done
    local d i; d="$(md5sum "$dir/checkpoints/horizyn_v1_0_dev.ckpt" | cut -d' ' -f1)"; i="$(md5sum "$dir/checkpoints/horizyn_v1_0_inf.ckpt" | cut -d' ' -f1)"
    [ "$d" = 5b1f938f8b0a82fbe91892a3b4e2bf2c ] && [ "$i" = cf6775b775287462099ae0681485a6bc ] || { echo "Horizyn checkpoint checksums differ (dev $d, inf $i)" >&2; exit 1; }
    record horizyn_dir "$dir" annotation; record horizyn_env "$he" annotation
    echo "Horizyn installed: $dir (conda environment '$he', checkpoints verified)"
}
install_colabfold() {
    need curl; need wget
    local dir="$TOOLS/localcolabfold"
    mkdir -p "$TOOLS" && cd "$TOOLS" || exit 1
    # localcolabfold's own installer (Miniforge + ColabFold 1.6 + JAX CUDA 12 + the AlphaFold2 weights); it installs into ./localcolabfold of the folder it runs in
    # (needs an NVIDIA driver only when ColabFold runs, not to install)
    [ -x "$dir/colabfold-conda/bin/colabfold_batch" ] || { fetch https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh install_colabbatch_linux.sh && bash install_colabbatch_linux.sh || exit 1; }
    cd "$REPO" || exit 1
    local cb="$dir/colabfold-conda/bin/colabfold_batch"
    "$cb" --help > /dev/null 2>&1 || { echo "colabfold_batch does not start" >&2; exit 1; }
    record colabfold_batch "$cb"
    echo "ColabFold installed: $cb"
}
install_chimerax() {
    need curl; need conda
    local pkg="" dir="$TOOLS/chimerax"
    while [ $# -gt 0 ]; do case "$1" in --chimerax-package) pkg="$2"; shift 2 ;; *) shift ;; esac; done
    [ -n "$pkg" ] || { echo "ChimeraX needs its licence accepted by a person: accept it at https://www.cgl.ucsf.edu/chimerax/download.html and pass the file or link you get:" >&2
                       echo "  bash workflow/setup_tools.sh --install chimerax --chimerax-package <ChimeraX-*.deb | .tar.gz | https://...>" >&2; exit 2; }
    mkdir -p "$dir/pkg" && cd "$dir/pkg" || exit 1
    case "$pkg" in http*) fetch "$pkg" chimerax_download || exit 1; pkg="$PWD/chimerax_download" ;; *) pkg="$(readlink -f "$pkg")" ;; esac
    case "$(file -b "$pkg")" in
        *Debian*|*debian*) dpkg-deb -x "$pkg" "$dir/root" || exit 1; bin="$(ls "$dir"/root/usr/lib/ucsf-chimerax*/bin/ChimeraX | head -n 1)" ;;
        *gzip*|*tar*)      mkdir -p "$dir/root" && tar xzf "$pkg" -C "$dir/root" || exit 1; bin="$(find "$dir/root" -path '*bin/ChimeraX' | head -n 1)" ;;
        *) echo "unrecognised ChimeraX package: $(file -b "$pkg")" >&2; exit 1 ;;
    esac
    [ -x "$bin" ] || { echo "no ChimeraX executable found in the package" >&2; exit 1; }
    # headless rendering needs libGLU.so.1 and libOpenGL.so.0; take them from conda-forge into a private folder instead of asking for root
    [ -d "$dir/libs" ] || conda create -y -p "$dir/libs" -c conda-forge libglu libglvnd > /dev/null || exit 1
    cat > "$dir/chimerax" << EOF
#!/usr/bin/env bash
export LD_LIBRARY_PATH="$dir/libs/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
exec "$bin" "\$@"
EOF
    chmod +x "$dir/chimerax"; cd "$REPO" || exit 1
    "$dir/chimerax" --offscreen --nogui --exit --cmd "version" || { echo "ChimeraX does not start headless" >&2; exit 1; }
    record chimerax "$dir/chimerax"
    echo "ChimeraX installed: $dir/chimerax"
}

case "$1" in
    --adopt) shift; for kv in "$@"; do k="${kv%%=*}"; v="${kv#*=}"; sec=tools; case "$k" in clean_*|horizyn_*) sec=annotation ;; esac; record "$k" "$v" "$sec"; done; exit 0 ;;
    --install)
        t="$2"; shift 2
        case "$t" in
            all) for x in foldseek clean horizyn colabfold; do (install_$x) || { echo "install of $x failed" >&2; exit 1; }; done ;;
            foldseek|clean|horizyn|colabfold) "install_$t" ;;
            chimerax) install_chimerax "$@" ;;
            *) echo "install: all | foldseek | clean | horizyn | colabfold | chimerax --chimerax-package <file|url>" >&2; exit 2 ;;
        esac
        exit $? ;;
esac

echo "Tools of the pipeline (refresh_config.local.json, refresh_config.json, tools/; EXTERNAL_TOOLS.md)"
for t in fimo streme mmseqs hmmbuild mafft; do command -v $t > /dev/null && say ok "$t (bioconda, from environment.yml)" || say miss "$t" "make env   (and conda activate smstpred)"; done
fs="$(get foldseek_bin)"; [ -z "$fs" ] && fs="$(command -v foldseek || true)"
[ -n "$fs" ] && say ok "foldseek: $fs" || say miss "foldseek (bake-off baseline, SaProt embeddings)" "bash workflow/setup_tools.sh --install foldseek"
cb="$(get colabfold_batch)"
[ -n "$cb" ] && say ok "colabfold_batch: $cb" || say miss "colabfold_batch (structures AlphaFold DB cannot supply)" "bash workflow/setup_tools.sh --install colabfold"
cx="$(get chimerax)"
[ -n "$cx" ] && say ok "ChimeraX: $cx" || say miss "ChimeraX (Figs 9 and 10 panels)" "accept the licence at https://www.cgl.ucsf.edu/chimerax/download.html, then: bash workflow/setup_tools.sh --install chimerax --chimerax-package <file or link>"
ca="$(get clean_app_dir)"
{ [ -n "$ca" ] && [ -f "$ca/CLEAN_infer_fasta.py" ]; } && say ok "CLEAN checkout: $ca" || say miss "CLEAN checkout (Stage 6)" "bash workflow/setup_tools.sh --install clean"
have_env "$(get clean_env)" && say ok "conda environment for CLEAN" || say miss "conda environment '$(get clean_env)'" "bash workflow/setup_tools.sh --install clean"
hz="$(get horizyn_dir)"
{ [ -n "$hz" ] && [ -s "$hz/checkpoints/horizyn_v1_0_inf.ckpt" ]; } && say ok "Horizyn checkout with checkpoints: $hz" || say miss "Horizyn checkout / checkpoints (Stage 6)" "bash workflow/setup_tools.sh --install horizyn"
have_env "$(get horizyn_env)" && say ok "conda environment for Horizyn" || say miss "conda environment '$(get horizyn_env)'" "bash workflow/setup_tools.sh --install horizyn"
echo; echo "$ok ok, $miss missing"
[ "$miss" -eq 0 ]
