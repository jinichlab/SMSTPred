#!/usr/bin/env python3
"""Run the registered steps (workflow/steps.py) on the release assembled from the archived record.

  python workflow/run_steps.py --mode verify --tier A      # CPU: every stage from the archived data, compared with the archived results
  python workflow/run_steps.py --mode verify --tier B      # GPU: re-train and compare outcomes with stated tolerances
  python workflow/run_steps.py --mode verify --only s2_fig2,s4_fig4
  python workflow/run_steps.py --mode run --stage 4        # just run, in place on the release folder (what `make stage4` does)

verify mode never touches the archived release: it copies the repository to a scratch folder and the release folder (symbolic links resolved) next to it, runs each
step there with RELEASE_DIR pointing at the copy, and compares every file the step wrote or changed with the archived file of the same relative path. The result
(PASS / FAIL / BUILT, the files compared, the log) goes to verify/report.json. The release folder is the newest release_* of $UNIPROT_PULL_DIR (default data/uniprot_pull),
made by `make data`; $RELEASE_DIR overrides it.
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
EXCLUDE_COPY = [".git", "data", "verify", "__pycache__", "*.pyc", ".verify_work", ".verify_release"]


def release_folder():
    sys.path.insert(0, str(REPO / "1_clustering"))
    from mmseqs_clusters import newest_release
    if os.environ.get("RELEASE_DIR"):
        return Path(os.environ["RELEASE_DIR"])
    return Path(newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))


def snapshot(root):
    snap = {}
    for p in root.rglob("*"):
        if p.is_file() and not p.is_symlink() and "__pycache__" not in p.parts:
            st = p.stat()
            snap[p.relative_to(root)] = (st.st_size, st.st_mtime_ns)
    return snap


def make_copies(work, ref_release, work_release):
    for d in (work, work_release):
        if d.exists():
            shutil.rmtree(d)
    subprocess.run(["rsync", "-a"] + [f"--exclude={e}" for e in EXCLUDE_COPY] + [f"{REPO}/", f"{work}/"], check=True)
    os.symlink(REPO / "data", work / "data")
    subprocess.run(["rsync", "-aL", f"{ref_release}/", f"{work_release}/"], check=True)      # -L: the archived big files are links into data/; the copy must be writable


def select(args):
    steps = STEPS
    if args.only:
        wanted = set(args.only.split(","))
        steps = [s for s in steps if s["id"] in wanted]
    tier = args.tier or (None if args.only else "A")
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
    args = ap.parse_args()

    steps = select(args)
    if args.list:
        for s in steps:
            print(f"{s['id']:<22s} tier {s['tier']}  stage {s['stage']:<3s} {s['cmd']}")
        return 0
    if not steps:
        print("no steps selected")
        return 1
    ref_release = release_folder()
    if not ref_release.exists():
        print(f"release folder {ref_release} not found: run `make data DATA_SRC=...` first")
        return 2

    verify = args.mode == "verify"
    work = Path(args.work) if verify else REPO
    work_release = Path(str(work) + "_release") if verify else ref_release
    logdir = REPO / "verify" / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    if verify:
        print(f"copying the repository to {work} and the release {ref_release.name} to {work_release} ...", flush=True)
        make_copies(work, ref_release, work_release)

    env = os.environ.copy()
    env.update(PYTHONHASHSEED="0", MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1", RELEASE_DIR=str(work_release), SMSTPRED_VERIFY="1")
    env.pop("RELEASE", None)
    results = []
    for s in steps:
        t0 = time.time()
        before = snapshot(work_release) if verify else {}
        log = logdir / f"{s['id']}.log"
        print(f"[{s['id']}] {s['cmd']}", flush=True)
        status, detail = "PASS", ""
        try:
            with open(log, "w") as fh:
                r = subprocess.run(["bash", "-c", s["cmd"]], cwd=work / s["cwd"], env=env, stdout=fh, stderr=subprocess.STDOUT, timeout=s["timeout"])
            if r.returncode != 0:
                status, detail = "FAIL", f"exit code {r.returncode} (see verify/logs/{s['id']}.log)"
        except subprocess.TimeoutExpired:
            status, detail = "FAIL", f"timed out after {s['timeout']} s"
        files = []
        if verify and status == "PASS" and s["check"]:
            files = [dict(path=what, status="numeric-match" if ok else "DIFFERENT", detail=d) for what, ok, d in CHECKS[s["check"]](work_release, ref_release)]
            bad = [f for f in files if f["status"] == "DIFFERENT"]
            if bad:
                status, detail = "FAIL", f"{len(bad)} of {len(files)} checks failed"
        elif verify and status == "PASS":
            after = snapshot(work_release)
            changed = sorted(p for p in after if before.get(p) != after[p] and p.suffix != ".log" and not any(str(p).endswith(x) for x in s["skip_compare"]))
            for p in changed:
                new, ref = work_release / p, ref_release / p
                if not ref.exists():
                    files.append(dict(path=str(p), status="new-file", detail="not in the archived record"))
                    continue
                st, d = compare(new, ref, rtol=s["rtol"], atol=s.get("atol", 1e-9), png_tol=s["png_tol"])
                files.append(dict(path=str(p), status=st, detail=d))
            bad = [f for f in files if f["status"] == "DIFFERENT"]
            if bad:
                status, detail = "FAIL", f"{len(bad)} file(s) differ from the archived version: " + ", ".join(f["path"] for f in bad[:4])
            else:
                matched = [f for f in files if f["status"] in ("identical", "numeric-match", "image-match")]
                if not matched:
                    status = "BUILT"
                    detail = "output written; no archived counterpart to compare" if files else "ran; wrote no files"
        dt = time.time() - t0
        results.append(dict(id=s["id"], stage=s["stage"], tier=s["tier"], paper=s["paper"], cmd=s["cmd"], status=status, detail=detail, seconds=round(dt, 1), files=files))
        print(f"    -> {status} ({dt:.0f} s) {detail}", flush=True)

    if verify:
        report = dict(generated=datetime.datetime.now().isoformat(timespec="seconds"), release=ref_release.name, versions=versions(), steps=results)
        out = REPO / "verify" / "report.json"
        prev = json.loads(out.read_text()) if out.exists() else None
        if prev and prev.get("release") == ref_release.name:      # a partial run updates only the steps it ran
            keep = {r["id"]: r for r in prev["steps"]}
            keep.update({r["id"]: r for r in results})
            report["steps"] = [keep[s["id"]] for s in STEPS if s["id"] in keep]
        out.write_text(json.dumps(report, indent=1))
        print(f"\nreport: {out}")
    n_fail = sum(r["status"] == "FAIL" for r in results)
    print(f"{sum(r['status'] == 'PASS' for r in results)} matched, {sum(r['status'] == 'BUILT' for r in results)} built without a reference, {n_fail} failed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
