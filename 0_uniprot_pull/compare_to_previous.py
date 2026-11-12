#!/usr/bin/env python3
"""
Stage 0 -- what changed between the archived data and a fresh pull.

  python compare_to_previous.py [--pull-dir DIR] [--release 2026_03] [--diagnose-removed 300]

Compares the newest release folder of the pull (pull_uniprot.py) with what the manuscript was built on:
  SAM set     new IPR029063 Swiss-Prot proteins vs the 10,619 proteins of embeddings_SAM_ProtT5.h5 (the proteins pulled for the original
              SAM set) and the 7,405 labelled ones in labels_SAM_full.csv; sequences that changed
  UP set      new reviewed enzymes vs labels_UP.csv (201,339); sequences that changed if the old FASTA is given (--old-up-fasta)
  proteomes   per organism: the new IPR029063 proteins (<= 1,000 aa) vs the proteins of 3_predictions/predictions_output/<organism>_*.csv
              (the 112,019 predicted ones), and the new proteome sizes / SAM counts vs figures/Fig5_proteome_scaling/proteome_vs_sam_bylabel.csv
  --diagnose-removed N   asks the API what happened to up to N removed SAM proteins (deleted, demoted to TrEMBL, lost the domain, longer ...)

Reads the archived tables from $SMSTPRED_DATA (or <repo>/data); writes comparison_report.txt and proteome_comparison.csv into the release folder.
"""
import argparse
import os
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

from uniprot_api import BASE, Client

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = Path(os.environ.get("SMSTPRED_DATA", REPO / "data"))
SAM = "IPR029063"
LINES = []


def say(msg=""):
    print(msg, flush=True)
    LINES.append(msg)


def newest_release(pull_dir):
    rel = sorted(p for p in Path(pull_dir).glob("release_*") if p.is_dir())
    if not rel:
        sys.exit(f"no release_* folder in {pull_dir}")
    return rel[-1]


def old_h5_keys(path):
    import h5py
    with h5py.File(path, "r") as f:
        return set(f.keys())


def set_report(name, old, new_all, new_le, old_label):
    kept = old & new_le
    gone = old - new_le
    gone_absent = gone - new_all
    say(f"  {old_label}: {len(old)} | new (all lengths): {len(new_all)}, <= length view: {len(new_le)}")
    say(f"    still present (<= length view): {len(kept)} ({len(kept) / len(old):.1%})")
    say(f"    removed: {len(gone)}: no longer in the new set at all {len(gone_absent)}, now longer than the length view {len(gone - gone_absent)}")
    say(f"    added (in the new set, not in the old one): {len(new_le - old)}")
    return gone


def diagnose(client, accessions, n):
    acc = sorted(accessions)[:n]
    found, cat = {}, Counter()
    for i in range(0, len(acc), 50):
        chunk = acc[i:i + 50]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        r = client.get(f"{BASE}/uniprotkb/search", {"query": q, "size": 50, "format": "tsv", "fields": "accession,reviewed,length,xref_interpro"})
        for line in r.text.splitlines()[1:]:
            a, rev, ln, ipr = (line.split("\t") + ["", "", "", ""])[:4]
            found[a] = (rev, int(ln or 0), SAM in ipr)
    for a in acc:
        if a not in found:
            cat["not found (deleted, merged or renamed)"] += 1
            continue
        rev, ln, has = found[a]
        if rev != "reviewed":
            cat["demoted to unreviewed"] += 1
        elif not has:
            cat["still reviewed, no longer annotated with IPR029063"] += 1
        elif ln > 1000:
            cat["reviewed, IPR029063, but longer than 1,000 aa"] += 1
        else:
            cat["reviewed, IPR029063, <= 1,000 aa (should be in the new set)"] += 1
    say(f"  what happened to {len(acc)} removed SAM proteins:")
    for k, v in cat.most_common():
        say(f"    {v:5d}  {k}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pull-dir")
    ap.add_argument("--release")
    ap.add_argument("--old-up-fasta", default=os.environ.get("UP_OLD_FASTA", ""), help="previous UniProt Swiss-Prot FASTA (or set UP_OLD_FASTA); optional")
    ap.add_argument("--diagnose-removed", type=int, default=0)
    args = ap.parse_args()
    pull_dir = Path(args.pull_dir or os.environ.get("UNIPROT_PULL_DIR") or REPO / "data" / "uniprot_pull")
    out = pull_dir / f"release_{args.release}" if args.release else newest_release(pull_dir)
    say(f"comparison of {out.name} with the archived data in {DATA}\n")
    L = int(__import__("json").loads((out / "pull_manifest.json").read_text()).get("length_view", 1000))

    # ---------------------------------------------------------------- SAM set
    say("== SAM set (Swiss-Prot, IPR029063)")
    new = pd.read_csv(out / "sam_mtases.tsv.gz", sep="\t", usecols=["Entry", "Length", "Sequence"], dtype=str)
    new["Length"] = new["Length"].astype(int)
    new_all, new_le = set(new.Entry), set(new.Entry[new.Length <= L])
    lab = pd.read_csv(DATA / "labels_SAM_full.csv", usecols=["Entry", "Sequence"], low_memory=False)
    gone = set()
    h5 = DATA / "embeddings_SAM_ProtT5.h5"
    if h5.exists():
        gone = set_report("pulled", old_h5_keys(h5), new_all, new_le, "the 10,619 proteins of the archived SAM embeddings")
    gone_lab = set_report("labelled", set(lab.Entry), new_all, new_le, "the 7,405 labelled training proteins")
    seq = dict(zip(new.Entry, new.Sequence))
    both = [(e, s) for e, s in zip(lab.Entry, lab.Sequence) if e in seq]
    changed = [e for e, s in both if str(s).strip().upper() != seq[e].strip().upper()]
    say(f"    sequence changed for {len(changed)} of the {len(both)} labelled proteins still present")
    if args.diagnose_removed and (gone or gone_lab):
        diagnose(Client(), gone or gone_lab, args.diagnose_removed)

    # ---------------------------------------------------------------- UP set
    say("\n== UP set (Swiss-Prot enzymes with catalytic activity)")
    up = pd.read_csv(out / "up_enzymes.tsv.gz", sep="\t", usecols=["Entry", "Length"] + (["Sequence"] if args.old_up_fasta and Path(args.old_up_fasta).exists() else []),
                     dtype=str)
    up["Length"] = up["Length"].astype(int)
    old_up = set(pd.read_csv(DATA / "labels_UP.csv", usecols=["Entry"]).Entry)
    set_report("", old_up, set(up.Entry), set(up.Entry[up.Length <= L]), "the 201,339 proteins of labels_UP.csv")
    if "Sequence" in up.columns:
        old_seq, cur = {}, None
        for line in open(args.old_up_fasta):
            if line.startswith(">"):
                cur = line[1:].split()[0]
                old_seq[cur] = []
            else:
                old_seq[cur].append(line.strip())
        nseq = dict(zip(up.Entry, up.Sequence))
        sh = [e for e in old_seq if e in nseq]
        ch = sum("".join(old_seq[e]) != nseq[e] for e in sh)
        say(f"    sequence changed for {ch} of the {len(sh)} proteins present in both ({args.old_up_fasta})")

    # ---------------------------------------------------------------- proteomes
    say("\n== proteomes (1,259 reference proteomes)")
    summ = pd.read_csv(out / "proteome_summary.csv")
    ok = summ[summ.status == "ok"]
    say(f"  ok {len(ok)}, empty {int((summ.status == 'empty').sum())}, errors {int(summ.status.astype(str).str.startswith('error').sum())} "
        f"of {len(summ)}; proteins {int(ok.n_all.sum())} (<= {L} aa {int(ok.n_le.sum())}); IPR029063 <= {L} aa {int(ok.n_ipr_le.sum())}")
    pred_dir = REPO / "3_predictions" / "predictions_output"
    rows = []
    for _, r in ok.iterrows():
        f = out / "reference_proteomes" / f"{r.stem}.csv"
        d = pd.read_csv(f, usecols=["Entry", "InterPro"], low_memory=False)
        new_set = set(d.Entry[d.InterPro.fillna("").str.contains(SAM, regex=False)])
        pf = pred_dir / f"{r.stem}_embeddings_predictions.csv"
        old_set = set(pd.read_csv(pf, usecols=["Entry"]).Entry) if pf.exists() else None
        rows.append({"organism": r.stem, "proteome_id": r.proteome_id, "n_proteins_new": r.n_le, "n_sam_new": len(new_set),
                     "n_sam_old_predicted": None if old_set is None else len(old_set),
                     "shared": None if old_set is None else len(new_set & old_set),
                     "added": None if old_set is None else len(new_set - old_set),
                     "removed": None if old_set is None else len(old_set - new_set)})
    cmp_df = pd.DataFrame(rows)
    cmp_df.to_csv(out / "proteome_comparison.csv", index=False)
    c = cmp_df.dropna(subset=["n_sam_old_predicted"])
    say(f"  organisms with an archived prediction file: {len(c)}")
    say(f"  IPR029063 proteins <= {L} aa: archived predicted set {int(c.n_sam_old_predicted.sum())} -> new {int(c.n_sam_new.sum())}; "
        f"shared {int(c.shared.sum())}, added {int(c.added.sum())}, removed {int(c.removed.sum())}")
    fig5 = REPO / "figures" / "Fig5_proteome_scaling" / "proteome_vs_sam_bylabel.csv"
    if fig5.exists():
        f5 = pd.read_csv(fig5).set_index("organism")
        j = cmp_df.set_index("organism").join(f5[["proteome_size", "n_sam"]], how="inner")
        say(f"  Fig 5 table (proteome_size, n_sam) vs new ({len(j)} organisms): proteome size {int(j.proteome_size.sum())} -> "
            f"{int(j.n_proteins_new.sum())} ({(j.n_proteins_new.sum() / j.proteome_size.sum() - 1):+.1%}); "
            f"n_sam {int(j.n_sam.sum())} -> {int(j.n_sam_new.sum())} ({(j.n_sam_new.sum() / j.n_sam.sum() - 1):+.1%})")
    (out / "comparison_report.txt").write_text("\n".join(LINES) + "\n")
    say(f"\nwrote {out / 'comparison_report.txt'} and proteome_comparison.csv")


if __name__ == "__main__":
    main()
