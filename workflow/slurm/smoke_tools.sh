#!/bin/bash
#SBATCH --job-name=smoke_tools
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=smoke_tools_%j.log

# Smoke test of the three external tools on a few proteins, through the same drivers the refresh chain uses, then a comparison with the archived results for those proteins:
#   CLEAN (EC numbers of 12 dark proteins), Horizyn (reaction scores of the same 12), ColabFold (one 70-residue sequence; needs api.colabfold.com).
#   REFERENCE=<release folder with dark_annotation/ clean/ and horizyn_run/ results>  sbatch workflow/slurm/smoke_tools.sh      (tools: workflow/setup_tools.sh)
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
[ -f Makefile ] || { echo "submit from the repository root" >&2; exit 1; }
eval "$(conda shell.bash hook)"; conda activate "${SMSTPRED_ENV:-smstpred}" || exit 1
REF="${REFERENCE:?set REFERENCE to a release folder holding the archived dark_annotation results}"
S="${SMOKE_DIR:-$PWD/.smoke_release}"; rm -rf "$S"; mkdir -p "$S/dark_annotation"
python - "$REF" "$S" << 'PY'
import sys, pandas as pd
from pathlib import Path
ref, s = Path(sys.argv[1]) / "dark_annotation", Path(sys.argv[2]) / "dark_annotation"
name = (ref / "dark_name.txt").read_text().strip()
m = pd.read_csv(ref / f"{name}_members.csv")
pick = pd.concat([m[m.cluster == c].head(2) for c in (317, 353, 512, 729, 517, 166)])
keep = set(pick.Entry)
pick.to_csv(s / f"{name}_members.csv", index=False)
(s / "dark_name.txt").write_text(name + "\n")
out, w = [], False
for line in open(ref / f"{name}.fasta"):
    if line.startswith(">"):
        w = line[1:].split()[0] in keep
    if w:
        out.append(line)
(s / f"{name}.fasta").write_text("".join(out))
print(len(keep), "proteins,", sum(l.startswith(">") for l in out), "sequences")
PY
export RELEASE_DIR="$S"
rc=0
echo "== CLEAN";   bash 6_dark_cluster_annotation/run_clean_release.sh   || { echo "CLEAN FAILED"; rc=1; }
echo "== Horizyn"; bash 6_dark_cluster_annotation/run_horizyn_release.sh || { echo "Horizyn FAILED"; rc=1; }
echo "== ColabFold"
CB="$(python workflow/cfgget.py colabfold_batch)"; mkdir -p "$S/cf"
printf ">smoke\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n" > "$S/cf/smoke.fasta"
"$CB" --num-models 1 --num-recycle 1 "$S/cf/smoke.fasta" "$S/cf/out" > "$S/cf/log.txt" 2>&1 && ls "$S/cf/out"/*.pdb > /dev/null 2>&1 && echo "ColabFold ok: $(ls "$S/cf/out"/*.pdb | head -1)" || { echo "ColabFold FAILED (see $S/cf/log.txt)"; tail -5 "$S/cf/log.txt"; rc=1; }
echo "== comparison with the archived results"
python - "$REF" "$S" << 'PY'
import sys, pandas as pd, numpy as np
from pathlib import Path
ref, s = Path(sys.argv[1]) / "dark_annotation", Path(sys.argv[2]) / "dark_annotation"
name = (s / "dark_name.txt").read_text().strip()
def rd(p):      # one line per protein: accession, then one or more EC calls with their distance (EC:x.y.z.w/0.1234)
    out = {}
    for line in open(p):
        f = line.strip().split(",")
        if len(f) >= 2:
            out[f[0]] = f[1].split("/")[0]
    return pd.Series(out)
a, b = rd(ref / "clean" / f"{name}_maxsep.csv"), rd(s / "clean" / f"{name}_maxsep.csv")
common = b.index.intersection(a.index)
print(f"CLEAN: {int((a[common] == b[common]).sum())} of {len(common)} proteins get the same top EC call as the archived run")
ta, tb = pd.read_csv(ref / "horizyn_run" / "horizyn_full_topk.csv"), pd.read_csv(s / "horizyn_run" / "horizyn_full_topk.csv")
x, y = ta[ta["rank"] == 1].set_index("accession").rhea_lr_id, tb[tb["rank"] == 1].set_index("accession").rhea_lr_id
c = y.index.intersection(x.index)
print(f"Horizyn: {int((x[c] == y[c]).sum())} of {len(c)} proteins have the same top-1 reaction as the archived run")
PY
exit $rc
