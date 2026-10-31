#!/usr/bin/env python3
"""Render verify/report.json as the pass/fail table in REPRODUCING.md (between the VERIFY-TABLE markers) and verify/report.md."""
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MARK = {"A": ("<!-- VERIFY-TABLE-START -->", "<!-- VERIFY-TABLE-END -->"), "B": ("<!-- VERIFY-GPU-TABLE-START -->", "<!-- VERIFY-GPU-TABLE-END -->")}
KIND = {"A": "files compared", "B": "checks"}


def table(report, tier="A"):
    steps = [s for s in report["steps"] if s["tier"] == tier]
    if not steps:
        return f"Run `make verify{'-gpu' if tier == 'B' else ''}` to fill this table.\n"
    rows = [f"| Step | Reproduces | Result | {KIND[tier].capitalize()} | Time |", "|---|---|---|---|---|"]
    for s in steps:
        counts = {}
        for f in s["files"]:
            counts[f["status"]] = counts.get(f["status"], 0) + 1
        if tier == "B":
            summary = f"{sum(f['status'] != 'DIFFERENT' for f in s['files'])} of {len(s['files'])} passed" if s["files"] else s.get("detail") or "-"
        else:
            summary = ", ".join(f"{n} {k}" for k, n in sorted(counts.items())) or s.get("detail") or "-"
        mark = {"PASS": "PASS", "BUILT": "built (no reference)"}.get(s["status"], "**FAIL**")
        rows.append(f"| `{s['id']}` | {s['paper']} | {mark} | {summary} | {s['seconds']:.0f} s |")
    if tier == "B":
        for s in steps:
            for f in s["files"]:
                cell = lambda x: str(x).replace("|", "\\|")
                rows.append(f"| &nbsp;&nbsp;`{s['id']}` check | {cell(f['path'])} | {'PASS' if f['status'] != 'DIFFERENT' else '**FAIL**'} | {cell(f['detail'])} | |")
    v = report["versions"]
    head = (f"Last verification: {report['generated']} (Python {v.get('python')}, numpy {v.get('numpy')}, pandas {v.get('pandas')}, "
            f"scipy {v.get('scipy')}, torch {v.get('torch')}, matplotlib {v.get('matplotlib')}). "
            f"{sum(s['status'] == 'PASS' for s in steps)} of {len(steps)} steps reproduce a committed result, "
            f"{sum(s['status'] == 'BUILT' for s in steps)} build their output without a committed reference to compare, "
            f"{sum(s['status'] == 'FAIL' for s in steps)} fail.\n")
    return head + "\n" + "\n".join(rows) + "\n"


def main():
    report = json.loads((REPO / "verify" / "report.json").read_text())
    md = REPO / "REPRODUCING.md"
    text = md.read_text()
    out = []
    for tier, (start, end) in MARK.items():
        t = table(report, tier)
        out.append(t)
        if start in text and end in text:
            a, b = text.index(start), text.index(end)
            text = text[:a + len(start)] + "\n" + t + text[b:]
            print(f"updated the tier {tier} table in REPRODUCING.md")
        else:
            print(f"markers {start} / {end} not found in REPRODUCING.md")
    md.write_text(text)
    (REPO / "verify" / "report.md").write_text("\n".join(out))


if __name__ == "__main__":
    main()
