#!/usr/bin/env python3
"""
Verify the apo/AF-motif upgrades: does each candidate substrate-bound structure's
protein actually carry the motif we want to place on it?

For each upgrade target, fetch the UniProt sequence, run FIMO with that class's
STREME motifs, and report which target motifs match and at what (UniProt) positions.
A >=50% cluster-identity link does NOT guarantee the exact motif is conserved, so
this is the gate before rebuilding any structure figure.

Outputs -> upgrade_verify_output/  (query FASTAs, FIMO results, summary.txt).
"""
import subprocess
import urllib.request
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "upgrade_verify_output"
OUT.mkdir(exist_ok=True)

# acc : (name, pdb, class, expected motif numbers)
TARGETS = {
    "Q8TEA1": ("NSUN6",     "5WWS", "RNA",     ["1", "4", "5"]),
    "P0C0R7": ("FtsJ/RlmE", "1EJ0", "RNA",     ["8"]),
    "Q9WVG6": ("CARM1",     "7OKP", "Protein", ["5"]),
}


def fetch(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    return urllib.request.urlopen(url, timeout=30).read().decode()


def main():
    by_class = {}
    for acc, (name, pdb, cls, exp) in TARGETS.items():
        fa = fetch(acc)
        print(f"[fetch] {acc} {name}  ({len(fa.splitlines())-1} seq lines)", flush=True)
        by_class.setdefault(cls, []).append((acc, fa))
    for cls, seqs in by_class.items():
        (OUT / f"query_{cls}.fasta").write_text("".join(fa for _, fa in seqs))

    # FIMO per class with the class's STREME motifs (same thresh as the original scan)
    for cls in by_class:
        subprocess.run(
            ["fimo", "--oc", str(OUT / f"fimo_{cls}"), "--thresh", "1e-4",
             str(HERE / f"streme_{cls}" / "streme.txt"), str(OUT / f"query_{cls}.fasta")],
            check=True)

    # parse + report
    lines = ["Upgrade-structure motif verification", "=" * 50]
    for acc, (name, pdb, cls, exp) in TARGETS.items():
        hits = {}
        with open(OUT / f"fimo_{cls}" / "fimo.tsv") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if not row or row[0].startswith("#") or row[0] == "motif_id":
                    continue
                mid, seqname, start, stop = row[0], row[2], row[3], row[4]
                if acc in seqname:
                    hits.setdefault(mid.split("-")[0], []).append((start, stop, mid))
        lines.append(f"\n{name}  ({acc}, PDB {pdb}, {cls})   want motif(s): {exp}")
        for mnum in exp:
            if mnum in hits:
                for start, stop, mid in hits[mnum]:
                    lines.append(f"   ✓ motif {mnum}  UniProt {start}-{stop}   [{mid}]")
            else:
                lines.append(f"   ✗ motif {mnum}  NOT found (below FIMO thresh)")
        other = sorted(set(hits) - set(exp), key=lambda x: int(x))
        if other:
            lines.append(f"   (also matches motifs: {', '.join(other)})")
    report = "\n".join(lines)
    print("\n" + report, flush=True)
    (OUT / "summary.txt").write_text(report + "\n")
    print(f"\n[written] {OUT}/summary.txt", flush=True)


if __name__ == "__main__":
    main()
