#!/usr/bin/env python3
"""Where is a tool?  python workflow/cfgget.py <key>      keys: clean_app_dir clean_env horizyn_dir horizyn_env colabfold_batch chimerax foldseek_bin

Order: the environment variable of the same name in capitals ($CLEAN_APP_DIR, $COLABFOLD_BATCH, ...), then refresh_config.local.json (this machine's paths, written by
`bash workflow/setup_tools.sh --install ...`, not committed), then refresh_config.json ("annotation" / "tools"), then the place `setup_tools.sh --install` puts the tool
(<repo>/tools/...). Prints nothing when the tool is not found. A path that does not exist is skipped, so a stale path of another machine in refresh_config.json never wins."""
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULTS = {
    "clean_app_dir": "tools/CLEAN/app",
    "horizyn_dir": "tools/horizyn",
    "colabfold_batch": "tools/localcolabfold/colabfold-conda/bin/colabfold_batch",
    "chimerax": "tools/chimerax/chimerax",
    "foldseek_bin": "tools/foldseek/bin/foldseek",
}
PLAIN = {"clean_env": "clean", "horizyn_env": "horizyn"}      # names, not paths


def load(name):
    p = REPO / name
    return json.loads(p.read_text()) if p.exists() else {}


def get(key):
    v = os.environ.get(key.upper())
    if v:
        return v
    if key in PLAIN:
        for c in (load("refresh_config.local.json"), load("refresh_config.json")):
            v = c.get("annotation", {}).get(key) or c.get("tools", {}).get(key)
            if v:
                return v
        return PLAIN[key]
    for c in (load("refresh_config.local.json"), load("refresh_config.json")):
        v = c.get("annotation", {}).get(key) or c.get("tools", {}).get(key)
        if v and Path(v).exists():
            return v
    d = REPO / DEFAULTS.get(key, "-")
    return str(d) if d.exists() else ""


if __name__ == "__main__":
    print(get(sys.argv[1]))
