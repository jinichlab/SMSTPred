#!/usr/bin/env python3
"""Run the registered steps, either in place (`--mode run`) or in a scratch copy with comparison (`--mode verify`).

verify mode never touches the committed files: it copies the repository to a scratch folder (data/ is linked, not copied), runs
each step there, and compares every file the step regenerated with the committed file of the same name. The result of each step
(PASS / FAIL, the files compared, the log) is written to verify/report.json.

  python workflow/run_steps.py --mode verify --tier A                  # everything runnable on CPU from the archived data
  python workflow/run_steps.py --mode verify --tier B                  # the GPU re-runs (retrain, compare outcomes with stated tolerances)
  python workflow/run_steps.py --mode verify --stage 2                 # one stage
  python workflow/run_steps.py --mode verify --only s2_table1,f_fig5   # chosen steps
  python workflow/run_steps.py --mode run --stage 7                    # just run, in place (what `make stage7` does)
"""
import argparse
import datetime
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from compare import compare            # noqa: E402
from steps import STEPS                # noqa: E402
from gpu_checks import CHECKS          # noqa: E402

REPO = HERE.parent
EXCLUDE_COPY = [".git", "data", "verify", "__pycache__", "*.pyc", ".verify_work"]


def snapshot(root):
    snap = {}
    for p in root.rglob("*"):
        rel = p.relative_to(root)
        if p.is_file() and not p.is_symlink() and "__pycache__" not in rel.parts and rel.parts[0] not in {"data", ".git", "verify"}:
            st = p.stat()
            snap[rel] = (st.st_size, st.st_mtime_ns)
    return snap


def make_work_copy(work: Path):
    if work.exists():
        shutil.rmtree(work)
    args = ["rsync", "-a"] + [f"--exclude={e}" for e in EXCLUDE_COPY] + [f"{REPO}/", f"{work}/"]
    subprocess.run(args, check=True)
    os.symlink(REPO / "data", work / "data")


def select(args):
    steps = STEPS
    if args.only:
        wanted = set(args.only.split(","))
        steps = [s for s in steps if s["id"] in wanted]
    tier = args.tier or (None if args.only else "A")        # the GPU re-runs (tier B) only run when asked for: they overwrite models in place
    if tier:
        steps = [s for s in steps if s["tier"] == tier]
    if args.stage:
        steps = [s for s in steps if s["stage"] == args.stage]
    return steps


def versions():
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for mod in ("numpy", "pandas", "scipy", "sklearn", "torch", "matplotlib", "PIL"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:
            pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["run", "verify"], default="verify")
    ap.add_argument("--tier")
    ap.add_argument("--stage")
    ap.add_argument("--only")
    ap.add_argument("--work", default=os.environ.get("SMSTPRED_VERIFY_WORK", str(REPO / ".verify_work")))
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--keep-going", action="store_true", default=True)
    args = ap.parse_args()

    steps = select(args)
    if args.list:
        for s in steps:
            print(f"{s['id']:<24s} tier {s['tier']}  stage {s['stage']:<8s} {s['cwd']}: {s['cmd']}")
        return 0
    if not steps:
        print("no steps selected")
        return 1

    verify = args.mode == "verify"
    root = Path(args.work) if verify else REPO
    logdir = REPO / "verify" / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    if verify:
        print(f"copying the repository to {root} ...", flush=True)
        make_work_copy(root)

    env = os.environ.copy()
    env.update(PYTHONHASHSEED="0", MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1")
    results = []
    for s in steps:
        t0 = time.time()
        before = snapshot(root) if verify else {}
        log = logdir / f"{s['id']}.log"
        print(f"[{s['id']}] {s['cmd']}", flush=True)
        status, detail = "PASS", ""
        try:
            with open(log, "w") as fh:
                r = subprocess.run(["bash", "-c", s["cmd"]], cwd=root / s["cwd"], env=env, stdout=fh, stderr=subprocess.STDOUT, timeout=s["timeout"])
            if r.returncode != 0:
                status, detail = "FAIL", f"exit code {r.returncode} (see verify/logs/{s['id']}.log)"
        except subprocess.TimeoutExpired:
            status, detail = "FAIL", f"timed out after {s['timeout']} s"
        files = []
        if verify and status == "PASS" and s["check"]:
            files = [dict(path=what, status="numeric-match" if ok else "DIFFERENT", detail=detail) for what, ok, detail in CHECKS[s["check"]](root, REPO)]
            bad = [f for f in files if f["status"] == "DIFFERENT"]
            if bad:
                status, detail = "FAIL", f"{len(bad)} of {len(files)} checks failed"
        elif verify and status == "PASS":
            after = snapshot(root)
            changed = sorted(p for p in after if before.get(p) != after[p])
            pairs = [(root / p, REPO / p) for p in changed]
            pairs += [(root / n, REPO / c) for n, c in s["compare_to"]]
            for new, ref in pairs:
                if not new.exists():
                    files.append(dict(path=str(ref.relative_to(REPO)), status="DIFFERENT", detail="expected output was not produced"))
                    continue
                st, d = compare(new, ref, rtol=s["rtol"], png_tol=s["png_tol"])
                files.append(dict(path=str(ref.relative_to(REPO)), status=st, detail=d))
            bad = [f for f in files if f["status"] == "DIFFERENT"]
            if bad:
                status, detail = "FAIL", f"{len(bad)} file(s) differ from the committed version"
            else:
                matched = [f for f in files if f["status"] in ("identical", "numeric-match", "image-match")]
                if not matched:
                    status = "BUILT"
                    detail = "output written; no committed counterpart to compare" if files else "ran; wrote no files"
        dt = time.time() - t0
        results.append(dict(id=s["id"], stage=s["stage"], tier=s["tier"], paper=s["paper"], cmd=s["cmd"], cwd=s["cwd"], status=status,
                            detail=detail, seconds=round(dt, 1), files=files))
        print(f"    -> {status} ({dt:.0f} s) {detail}", flush=True)

    if verify:
        report = dict(generated=datetime.datetime.now().isoformat(timespec="seconds"), versions=versions(), steps=results)
        out = REPO / "verify" / "report.json"
        prev = json.loads(out.read_text()) if out.exists() else None
        if prev:   # a partial run updates only the steps it ran
            keep = {r["id"]: r for r in prev["steps"]}
            keep.update({r["id"]: r for r in results})
            order = [s["id"] for s in STEPS]
            report["steps"] = [keep[i] for i in order if i in keep]
        out.write_text(json.dumps(report, indent=1))
        print(f"\nreport: {out}")
    n_fail = sum(r["status"] == "FAIL" for r in results)
    print(f"{sum(r['status'] == 'PASS' for r in results)} matched, {sum(r['status'] == 'BUILT' for r in results)} built without a reference, {n_fail} failed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
