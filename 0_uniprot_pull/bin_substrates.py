#!/usr/bin/env python3
"""
Stage 0 -- substrate typing of a UniProt pull: Metabolite / Protein / RNA / DNA labels from the reaction annotation.

  python bin_substrates.py sam|up|all [--pull-dir DIR] [--release 2026_03] [--benchmark-max-len 1000] [--deployment-max-len 2000]
                                      [--benchmark-sources reaction_curated] [--deployment-sources reaction_curated]

Reads sam_mtases.tsv.gz (mode sam: only the SAM reactions are typed) and / or up_enzymes.tsv.gz (mode all: every reaction) of the newest
release folder written by pull_uniprot.py, applies the deterministic rules of substrate_binning.py and writes next to them:

  labels_SAM.csv / labels_UP.csv   one row per protein: the four 0/1 label columns (same names as the archived labels_SAM_full.csv / labels_UP.csv),
                                   label_source (reaction_curated | reaction_uncurated | name_rule | none), the SAM transfer kind, the Rhea ids,
                                   the acceptors with the rule that typed each one, and the two length flags:
                                     in_benchmark   Length <= --benchmark-max-len  (1,000) and label_source in --benchmark-sources
                                     in_deployment  Length <= --deployment-max-len (2,000) and label_source in --deployment-sources
  binning_rules_used_*.csv         every rule that fired (Rhea family / word rule; polypeptide polymers by name) with the number of proteins: the review list
  binning_text_fallback_names_*.csv  the participant names that Rhea's table did not type (no Rhea id in the UniProt reaction) and the word rules did
  binning_vs_archived.csv          the proteins whose label set differs from the archived labelling, with both sets and the acceptors
  binning_report.txt               counts, the agreement with the archived labels, and the measured accuracy of the weaker name-rule tier
  binning_manifest.json            inputs (sha256), parameters, rules version

Participants are typed from Rhea's own compound families (rhea_compounds.tsv.gz, written by fetch_rhea_compounds.py; the copy next to this script is the
default, --rhea-table overrides it, --refresh-rhea downloads a new one first). Polypeptide polymer names never need a rule: every polypeptide compound
Rhea has ever described is a protein-bound residue. The one list that is curated by name is the set of polymer names already reviewed
(polypeptide_tokens_reviewed.txt): a run prints the polypeptide polymers it has not seen before, so that a new redox partner or carrier
that Rhea added since is looked at once (--accept-new-tokens adds them).

Reads the archived labels from $SMSTPRED_DATA (default <repo>/data). Pure CPU, about a minute for the UP set.
"""
import argparse
import hashlib
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd

import substrate_binning as sb

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = Path(os.environ.get("SMSTPRED_DATA", REPO / "data"))
PULL = Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull"))
RULES_VERSION = "2026-11-02"
LINES = []


def say(msg=""):
    print(msg, flush=True)
    LINES.append(msg)


def newest_release(pull_dir, release=None):
    if release:
        return Path(pull_dir) / f"release_{release}"
    rel = sorted(p for p in Path(pull_dir).glob("release_*") if p.is_dir() and not p.name.startswith("_"))
    if not rel:
        sys.exit(f"no release_* folder in {pull_dir}")
    return rel[-1]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def label_set(row):
    return frozenset(c for c in sb.CLASSES if row[c] == 1)


def rule_key(cls, rule):
    """The rule without the matched word, so that every RNA word counts as one rule; bracket / polymer rules keep their token."""
    for w in ("RNA word", "DNA word"):
        if f"{w}:" in rule:  # a text rule 'text:RNA word:tRNA' -> 'text:RNA word'; the structural 'rhea:Polynucleotide (RNA word)' has no matched word to drop
            return cls, rule.split(f"{w}:")[0] + w
    return cls, rule


def bin_table(df, mode, rhea):
    """-> (df with the four label columns, label_source, transfer, n_reactions, RheaID, Acceptors; Counter of the rules that fired, one per protein;
    Counter of the names typed by the word rules, one per protein)."""
    rows, rules, fallback = [], Counter(), Counter()
    for ca, names, ec, fn in zip(df["Catalytic activity"], df["Protein names"], df["EC number"], df["Function [CC]"]):
        r = sb.bin_entry(ca, mode, rhea)
        labels, source = set(r["labels"]), r["source"]
        if source == "none" and mode == "sam":  # no typable reaction: the weaker tier, from the names (written for methyltransferase names; not used for UP)
            name_labels, _ = sb.label_from_names(names, ec, fn)
            if name_labels:
                labels, source = set(name_labels), "name_rule"
        rules.update({rule_key(c, rl) for _, c, rl in r["typed"]})
        fallback.update({(n[:90], c, rl) for n, c, rl in r["typed"] if rl.startswith(("text:", "rhea:Polymer>"))})
        rows.append({**{c: int(c in labels) for c in sb.CLASSES}, "label_source": source, "transfer": ";".join(sorted(r["transfer"])),
                     "Protein_via_carrier": int(r["carrier_only"]), "n_reactions": r["n_reactions"], "n_uncurated_ignored": r["n_ignored"], "RheaID": ";".join(r["rhea"]), "Acceptors": " | ".join(r["acceptors"])})
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1), rules, fallback


def rules_table(rules):
    """How many proteins each rule typed (one count per protein and rule); the review list of the rule set."""
    return pd.DataFrame([(c, r, n) for (c, r), n in rules.items()], columns=["class", "rule", "n_proteins"]).sort_values(["class", "n_proteins"], ascending=[True, False])


def rhea_table_audit(rhea):
    """Which of the two safeguards against Rhea's mis-filings the table still needs, so that a Rhea correction is noticed and not silently absorbed."""
    mis = sorted(n for n, fam in rhea.items() if fam == "Polypeptide" and (sb.RNA_RE.search(sb._outside_brackets(n)) or sb.DNA_RE.search(sb._outside_brackets(n))))
    say(f"  Rhea table check: nucleic-acid residues still filed as polypeptide (the conflict rule retypes them): {len(mis)}" + (" -> Rhea has corrected them all, the rule is dormant" if not mis else ""))
    for n in mis:
        say(f"      {n[:100]}")
    for n, fam in sb.FAMILY_OVERRIDES.items():
        cur = rhea.get(n)
        say(f"  Rhea table check: override {n!r} -> {fam}; Rhea now files it as {cur}" + (" (Rhea agrees: the override can be removed)" if cur == fam else "" if cur else " (not in the table any more: renamed or removed)"))


def structural_share(rules):
    """Share of the (protein, rule) typings that came from Rhea's families rather than from the word rules."""
    tot = sum(rules.values())
    rhea = sum(n for (_, r), n in rules.items() if r.startswith("rhea:") and not r.startswith("rhea:Polymer>"))
    poly = sum(n for (_, r), n in rules.items() if r.startswith("rhea:Polymer>"))
    text = sum(n for (_, r), n in rules.items() if r.startswith("text:"))
    say(f"  acceptor typings (one per protein and rule): {tot} | from Rhea's compound families {rhea} ({100 * rhea / max(tot, 1):.1f}%) | Rhea Polymer family, then word rules {poly} "
        f"| names Rhea does not know, word rules {text} ({100 * text / max(tot, 1):.1f}%)")


def new_polymer_tokens(rules, reviewed_path, accept):
    """The polypeptide polymers (bracket names) that typed a Protein and are not in the reviewed list; a new redox partner would show up here."""
    seen = {r.split(":", 1)[1] for (_, r), _ in rules.items() if r.startswith("rhea:Polypeptide:")}
    seen |= {r.split(":", 1)[1] for (_, r), _ in rules.items() if r.startswith("rhea:Polymer>") and ":[" in r}
    reviewed = set(Path(reviewed_path).read_text().splitlines()) if Path(reviewed_path).exists() else None
    if reviewed is None:
        Path(reviewed_path).write_text("\n".join(sorted(seen)) + "\n")
        say(f"  polymer names: {len(seen)} written as the reviewed list {Path(reviewed_path).name} (first run: review it)")
        return
    new = sorted(seen - reviewed)
    if not new:
        say(f"  polymer names: all {len(seen)} that typed a Protein are in the reviewed list")
        return
    say(f"  NEW polymer names that typed a Protein since the reviewed list ({len(new)}): look at these once (a redox partner or carrier belongs in REDOX_PARTNER):")
    for t in new:
        say(f"      {t}")
    if accept:
        Path(reviewed_path).write_text("\n".join(sorted(reviewed | set(new))) + "\n")
        say("  (added to the reviewed list: --accept-new-tokens)")


def flags(df, bench_max, dep_max, bench_sources, dep_sources):
    labelled = df[list(sb.CLASSES)].sum(axis=1) > 0
    df["in_benchmark"] = (df.Length <= bench_max) & df.label_source.isin(bench_sources) & labelled
    df["in_deployment"] = (df.Length <= dep_max) & df.label_source.isin(dep_sources) & labelled
    return df


def summarize(df, name, bench_max, dep_max):
    say(f"\n== {name}: {len(df)} proteins")
    say(f"  label_source: {dict(df.label_source.value_counts())}")
    for src in ("reaction_curated", "reaction_uncurated", "name_rule"):
        d = df[df.label_source == src]
        if len(d):
            say(f"    {src:19s} {len(d):7d} | " + ", ".join(f"{c} {int(d[c].sum())}" for c in sb.CLASSES) + f" | multi-label {int((d[list(sb.CLASSES)].sum(axis=1) > 1).sum())}")
    if "transfer" in df and df.transfer.ne("").any():
        say(f"  SAM transfer kind (reaction-typed): {dict(df[df.transfer != ''].transfer.value_counts())}")
    pc = df[(df.Protein == 1) & df.label_source.str.startswith("reaction")]
    say(f"  Protein-labelled enzymes: {len(pc)}, of which only through a carrier-tethered substrate (ACP, sulfur / lipoyl / peptidyl carrier ...; column Protein_via_carrier): {int(pc.Protein_via_carrier.sum())}")
    say(f"  in_benchmark (<= {bench_max} aa): {int(df.in_benchmark.sum())} | in_deployment (<= {dep_max} aa): {int(df.in_deployment.sum())}")
    for flag in ("in_benchmark", "in_deployment"):
        d = df[df[flag]]
        say(f"    {flag:14s} " + ", ".join(f"{c} {int(d[c].sum())}" for c in sb.CLASSES))


def vs_archived(df, archived_path, name, out_dir, key="Entry"):
    if not Path(archived_path).exists():
        say(f"\n== {name} vs archived: {archived_path} not found, skipped")
        return
    old = pd.read_csv(archived_path, low_memory=False, usecols=[key, *sb.CLASSES])
    cols = list(sb.CLASSES)
    both = df.merge(old, on=key, suffixes=("", "_old"))
    typed = both[both[cols].sum(axis=1) > 0]
    say(f"\n== {name} vs the archived labels ({archived_path.name}, {len(old)} proteins)")
    say(f"  archived proteins in the new pull: {len(both)} ({100 * len(both) / len(old):.1f}%) | typed by the new rules: {len(typed)} "
        f"(source: {dict(typed.label_source.value_counts())})")
    new_set = typed.apply(label_set, axis=1)
    old_set = typed.apply(lambda r: frozenset(c for c in sb.CLASSES if r[c + "_old"] == 1), axis=1)
    agree = new_set == old_set
    say(f"  identical label sets: {int(agree.sum())} of {len(typed)} ({100 * agree.mean():.2f}%)")
    for c in sb.CLASSES:
        o, n = old_set.map(lambda s: c in s), new_set.map(lambda s: c in s)
        say(f"    {c:11s} archived {int(o.sum()):6d} | new rules {int(n.sum()):6d} | both {int((o & n).sum()):6d} | archived only {int((o & ~n).sum()):4d} | new only {int((~o & n).sum()):4d}")
    dis = typed[~agree.values].copy()
    dis["archived_labels"] = [",".join(sorted(s)) for s in old_set[~agree]]
    dis["new_labels"] = [",".join(sorted(s)) for s in new_set[~agree]]
    pat = Counter(zip(dis.archived_labels, dis.new_labels)).most_common(8)
    say(f"  differing proteins by (archived -> new): {pat}")
    untyped = both[both[cols].sum(axis=1) == 0]
    say(f"  archived proteins the new rules leave untyped: {len(untyped)}")
    keep = [c for c in (key, "Protein names", "Length", "archived_labels", "new_labels", "label_source", "Acceptors") if c in dis]
    dis[keep].to_csv(out_dir / "binning_vs_archived.csv", index=False)
    new_only = df[(~df[key].isin(old[key])) & (df.label_source != "none")]
    say(f"  typed by the new rules, not in the archived set: {len(new_only)} (source: {dict(new_only.label_source.value_counts())}); "
        f"new reaction-typed ones <= 1,000 aa: {int(((new_only.label_source.str.startswith('reaction')) & (new_only.Length <= 1000)).sum())}")


def name_rule_accuracy(df):
    """Accuracy of label_from_names on proteins typed from a reaction (the only labels it can be scored against). Two references: the curated
    reactions (well-characterised proteins, informative names: flattering) and the by-similarity ones (less studied: the fairer estimate)."""
    for src in ("reaction_curated", "reaction_uncurated"):
        ref = df[df.label_source == src].reset_index(drop=True)
        if not len(ref):
            continue
        names = [sb.label_from_names(n, e, f)[0] for n, e, f in zip(ref["Protein names"], ref["EC number"], ref["Function [CC]"])]
        fired = [i for i, s in enumerate(names) if s]
        say(f"\n== accuracy of the name-only tier (the name rule never sees the reaction), scored on {len(ref)} {src} proteins")
        say(f"  fires on {len(fired)} ({100 * len(fired) / len(ref):.1f}%)")
        for c in ("RNA", "DNA", "Protein"):
            pred = [i for i in fired if c in names[i]]
            right = sum(1 for i in pred if ref.loc[i, c] == 1)
            truth = int(ref[c].sum())
            say(f"    predicts {c:8s}: {len(pred):6d} | precision {100 * right / max(len(pred), 1):5.1f}% | recall of all {c} {100 * right / max(truth, 1):5.1f}% ({right} of {truth})")
        meta = set(ref.index[(ref.Metabolite == 1) & (ref[["RNA", "DNA", "Protein"]].sum(axis=1) == 0)])
        wrong = sum(1 for i in fired if i in meta)
        say(f"    of {len(meta)} metabolite-only proteins the name rule wrongly fires on {wrong} ({100 * wrong / max(len(meta), 1):.1f}%)")


def run(which, args):
    rel = newest_release(args.pull_dir, args.release)
    out = Path(args.out) if args.out else rel
    out.mkdir(parents=True, exist_ok=True)
    say(f"substrate typing of {rel} (rules {RULES_VERSION}), output {out}")
    rhea_path = Path(args.rhea_table)
    if args.refresh_rhea:
        import fetch_rhea_compounds
        sys.argv = ["fetch_rhea_compounds.py", "--out", str(rhea_path)]
        fetch_rhea_compounds.main()
    if not rhea_path.exists():
        sys.exit(f"{rhea_path} not found: run fetch_rhea_compounds.py (or pass --refresh-rhea)")
    rhea = sb.load_rhea_table(rhea_path)
    say(f"Rhea compound table {rhea_path.name}: {len(rhea)} compounds, sha256 {sha256(rhea_path)[:16]}")
    rhea_table_audit(rhea)
    (out / "rhea_compounds.tsv.gz").write_bytes(rhea_path.read_bytes()) if rhea_path.resolve() != (out / "rhea_compounds.tsv.gz").resolve() else None
    manifest = {"rules_version": RULES_VERSION, "rhea_table_sha256": sha256(rhea_path), "release_dir": str(rel), "benchmark_max_len": args.benchmark_max_len, "deployment_max_len": args.deployment_max_len,
                "benchmark_sources": args.benchmark_sources, "deployment_sources": args.deployment_sources, "inputs": {}, "outputs": {}}
    parts = {"sam": ("sam_mtases.tsv.gz", "labels_SAM.csv", "sam", "SAM MTases (IPR029063)", DATA / "labels_SAM_full.csv"),
             "up": ("up_enzymes.tsv.gz", "labels_UP.csv", "all", "Swiss-Prot enzymes (UP)", DATA / "labels_UP.csv")}
    for w in (["sam", "up"] if which == "all" else [which]):
        fname, oname, mode, title, archived = parts[w]
        src = rel / fname
        if not src.exists():
            sys.exit(f"{src} not found: run pull_uniprot.py {w} first")
        t0 = time.time()
        use = ["Entry", "Protein names", "Length", "Annotation", "EC number", "Function [CC]", "Catalytic activity"]
        df = pd.read_csv(src, sep="\t", usecols=lambda c: c in use + ["Entry Name"], low_memory=False)
        df, rules, fallback = bin_table(df, mode, rhea)
        df = flags(df, args.benchmark_max_len, args.deployment_max_len, args.benchmark_sources, args.deployment_sources)
        keep = ["Entry", "Entry Name", "Protein names", "Length", "Annotation", "EC number", *sb.CLASSES, "label_source", "transfer", "Protein_via_carrier", "n_reactions", "n_uncurated_ignored", "RheaID", "Acceptors",
                "in_benchmark", "in_deployment"]
        df[[c for c in keep if c in df]].to_csv(out / oname, index=False)
        rules_table(rules).to_csv(out / f"binning_rules_used_{w}.csv", index=False)
        pd.DataFrame([(n, c, rl, k) for (n, c, rl), k in fallback.most_common()], columns=["name", "class", "rule", "n_proteins"]).to_csv(
            out / f"binning_text_fallback_names_{w}.csv", index=False)
        summarize(df, title, args.benchmark_max_len, args.deployment_max_len)
        structural_share(rules)
        new_polymer_tokens(rules, HERE / "polypeptide_tokens_reviewed.txt", args.accept_new_tokens)
        if mode == "sam":
            name_rule_accuracy(df)
        if not args.no_compare:
            vs_archived(df, archived, title, out)
            if (out / "binning_vs_archived.csv").exists():
                (out / "binning_vs_archived.csv").rename(out / f"binning_vs_archived_{w}.csv")
        manifest["inputs"][fname] = sha256(src)
        manifest["outputs"][oname] = {"rows": len(df), "seconds": round(time.time() - t0, 1)}
    (out / "binning_report.txt").write_text("\n".join(LINES) + "\n")
    (out / "binning_manifest.json").write_text(json.dumps(manifest, indent=2))
    say(f"\nwrote {out / 'binning_report.txt'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("which", choices=["sam", "up", "all"])
    ap.add_argument("--pull-dir", default=str(PULL))
    ap.add_argument("--release")
    ap.add_argument("--out", help="output folder (default: the release folder)")
    ap.add_argument("--benchmark-max-len", type=int, default=1000)
    ap.add_argument("--deployment-max-len", type=int, default=2000)
    ap.add_argument("--benchmark-sources", nargs="+", default=["reaction_curated"], choices=["reaction_curated", "reaction_uncurated", "name_rule"])
    ap.add_argument("--deployment-sources", nargs="+", default=["reaction_curated"], choices=["reaction_curated", "reaction_uncurated", "name_rule"])
    ap.add_argument("--no-compare", action="store_true", help="skip the comparison with the archived labels")
    ap.add_argument("--rhea-table", default=str(HERE / "rhea_compounds.tsv.gz"), help="Rhea compound families (default: the copy next to this script)")
    ap.add_argument("--refresh-rhea", action="store_true", help="download a new Rhea compound table first (network)")
    ap.add_argument("--accept-new-tokens", action="store_true", help="add the polymer names that are new since the reviewed list to it")
    args = ap.parse_args()
    run(args.which, args)


if __name__ == "__main__":
    main()
