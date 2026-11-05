#!/usr/bin/env python3
"""Compare a release folder rebuilt by the pipeline with an archived release folder, file by file.

  python workflow/compare_releases.py --new <rebuilt release> --ref <archived release> [--rtol 1e-3] [--out compare_report.csv]

Every file that exists in both folders is compared with the rules of `make verify` (workflow/compare.py): tables and statistics number by number within --rtol, figures by
pixel difference, binary models / arrays / archives not compared. Prints the count per status and lists the files that differ, then the key numbers side by side (key_numbers.md).
Tolerances are loose on purpose (default 1e-3): this compares two independent rebuilds on the same inputs, where retrained networks and library versions move numbers slightly."""
import argparse
import difflib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import compare            # noqa: E402

SKIP_DIRS = {"structures", "reference_proteomes", "reference_proteomes_long", "_smoke", "tmp_mmseqs", "colabfold_output", "cif", "predictions", "proteomes_ProtT5", "article_figures", "horizyn_run", "clean_run"}
SKIP_NAMES = {"context_cache.json", "proteome_status.jsonl", "colabfold_run.txt"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--rtol", type=float, default=1e-3)
    ap.add_argument("--out", default="compare_report.csv")
    a = ap.parse_args()
    new, ref = Path(a.new), Path(a.ref)
    rows = []
    for p in sorted(new.rglob("*")):
        if not p.is_file() or p.is_symlink():
            continue
        rel = p.relative_to(new)
        if set(rel.parts[:-1]) & SKIP_DIRS or rel.name in SKIP_NAMES or rel.name.endswith((".log", ".pyc")) or p.stat().st_size > 400e6:
            continue
        q = ref / rel
        if not q.exists():
            rows.append((str(rel), "only-in-new", ""))
            continue
        try:
            st, d = compare(p, q, rtol=a.rtol, png_tol=12.0)
        except Exception as e:      # an unreadable file is a finding, not a crash
            st, d = "DIFFERENT", f"could not compare: {type(e).__name__}"
        rows.append((str(rel), st, str(d)[:200]))
    df = pd.DataFrame(rows, columns=["file", "status", "detail"])
    df.to_csv(a.out, index=False)
    print(df.status.value_counts().to_string())
    bad = df[df.status == "DIFFERENT"]
    print(f"\n{len(bad)} files differ:")
    for r in bad.itertuples():
        print(f"  {r.file}: {r.detail}")
    kn, kr = new / "key_numbers.md", ref / "key_numbers.md"
    if kn.exists() and kr.exists():
        print("\nkey_numbers.md (- archived, + rebuilt):")
        for l in difflib.unified_diff(kr.read_text().splitlines(), kn.read_text().splitlines(), lineterm="", n=0):
            if not l.startswith(("---", "+++", "@@")):
                print(" ", l[:230])


if __name__ == "__main__":
    main()
