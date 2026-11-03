#!/usr/bin/env python3
"""Assemble the FIMO cross-scan (run_fimo_cross.sh) into interpretable tables.

The per-class diagonal (fimo_all_*) could only ever report a motif's coverage
inside its own class. Scanning every motif against every protein lets us ask the
question the "discriminative motif" claim actually rests on: how often does each
motif occur OUTSIDE the class it was discovered in?

Outputs (this dir):
  motif_specificity_matrix.csv   per motif: proteins carrying it in EVERY class,
                                 as counts and as class-size-normalised rates,
                                 plus TWO specificity ratios and a significance test.
                                 Rates matter because classes differ ~6x in size.
  fimo_cross_protein_motifs.csv  per protein: class, which motifs it carries, and
                                 how many are own-class vs foreign-class. Includes
                                 the Mixed/Unassigned proteins the diagonal skipped.
  fimo_cross_summary.txt         human-readable report.

TWO SPECIFICITY RATIOS (they disagree, and the disagreement is informative)
  ratio_vs_max  = own_pct / (highest single off-class rate).  Conservative against a
                  motif concentrated in ONE rival class, but blind to a motif that is
                  diffusely common everywhere (40% own vs 20/20/20 scores 2.0 despite
                  being in a fifth of everything).
  ratio_vs_bg   = own_pct / (rate across ALL non-own proteins pooled).  Penalises the
                  diffuse case correctly. This is the primary metric.

SIGNIFICANCE
  Fisher exact test, own-class carriers vs pooled-background carriers, Bonferroni-
  corrected across the 32 motifs, plus a 95% CI on the log2 odds ratio (Woolf).
  With n up to 118,816 nearly everything is "significant", so the CI/effect size --
  not the p-value -- is what should drive any cutoff. The p-value only guards against
  calling a small-n class (DNA, n=6,822) specific on noise.
"""
import argparse
import math
from pathlib import Path
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CROSS = HERE / "fimo_cross"
CLASSES = ["DNA", "RNA", "Protein", "Metabolite"]
N_MOTIFS = 32          # Bonferroni family size


def load_hits():
    """Concatenate the four cross-scan TSVs, tagging each motif with its origin class."""
    frames = []
    for cls in CLASSES:
        f = CROSS / f"{cls}.tsv"
        if not f.exists():
            print(f"  [warn] missing {f}")
            continue
        d = pd.read_csv(f, sep="\t", comment="#")
        d = d.dropna(subset=["sequence_name"])
        d["motif_class"] = cls
        # motif_id looks like "2-RJVYSTCS" -> number 2, consensus RJVYSTCS
        d["motif_num"] = d["motif_id"].astype(str).str.split("-").str[0]
        d["motif"] = cls + "-" + d["motif_num"]
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def main():
    global HERE, CROSS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=str(HERE), help="folder with corpus.csv and fimo_cross/ (a release: <release>/motif_analysis); outputs are written there")
    HERE = Path(ap.parse_args().dir).resolve()
    CROSS = HERE / "fimo_cross"
    corpus = pd.read_csv(HERE / "corpus.csv", usecols=["Entry", "class"],
                         keep_default_na=False, low_memory=False)
    corpus["class"] = corpus["class"].replace({"None": "Unassigned"})
    cls_of = dict(zip(corpus["Entry"], corpus["class"]))
    sizes = corpus["class"].value_counts().to_dict()
    all_classes = ["DNA", "RNA", "Protein", "Metabolite", "Mixed", "Unassigned"]

    hits = load_hits()
    hits["prot_class"] = hits["sequence_name"].map(cls_of)
    unknown = hits["prot_class"].isna().sum()
    if unknown:
        print(f"  [warn] {unknown} hits on sequences absent from corpus.csv")
        hits = hits.dropna(subset=["prot_class"])
    print(f"[cross] {len(hits):,} hits over {hits.sequence_name.nunique():,} proteins")

    # ---- 1. specificity matrix: distinct proteins per (motif, protein class) ----
    pc = (hits.drop_duplicates(["motif", "sequence_name"])
              .groupby(["motif", "prot_class"]).size().unstack(fill_value=0))
    for c in all_classes:
        if c not in pc.columns:
            pc[c] = 0
    pc = pc[all_classes]

    rows = []
    cons = (hits.drop_duplicates("motif").set_index("motif")["motif_id"]
                .astype(str).str.split("-").str[1].to_dict())
    for motif, r in pc.iterrows():
        own = motif.split("-")[0]
        rec = {"motif": motif, "consensus": cons.get(motif, ""), "own_class": own}
        for c in all_classes:
            rec[f"n_{c}"] = int(r[c])
            rec[f"pct_{c}"] = round(100 * r[c] / sizes.get(c, 1), 2)
        others = [c for c in CLASSES if c != own]
        rec["own_pct"] = rec[f"pct_{own}"]

        # --- ratio 1: vs the single worst off-class (conservative but blind to diffuse) ---
        rec["max_offclass_pct"] = max(rec[f"pct_{c}"] for c in others)
        rec["worst_offclass"] = max(others, key=lambda c: rec[f"pct_{c}"])
        rec["ratio_vs_max"] = (round(rec["own_pct"] / rec["max_offclass_pct"], 2)
                               if rec["max_offclass_pct"] > 0 else float("inf"))

        # --- ratio 2: vs pooled background = every protein NOT of the own class ---
        # (includes Mixed/Unassigned: they are real proteins that carry the motif or
        # not, and a motif common across them is not class-specific.)
        a = int(r[own])                                   # own-class carriers
        n_own = sizes.get(own, 0)
        b = n_own - a                                     # own-class non-carriers
        c_ = int(sum(r[x] for x in all_classes if x != own))   # background carriers
        n_bg = sum(sizes.get(x, 0) for x in all_classes if x != own)
        d = n_bg - c_                                     # background non-carriers
        rec["n_background"] = c_
        rec["background_pct"] = round(100 * c_ / n_bg, 2) if n_bg else 0.0
        rec["ratio_vs_bg"] = (round(rec["own_pct"] / rec["background_pct"], 2)
                              if rec["background_pct"] > 0 else float("inf"))

        # --- significance: Fisher exact (own vs pooled background), Bonferroni ---
        odds, p = stats.fisher_exact([[a, b], [c_, d]], alternative="two-sided")
        rec["odds_ratio"] = round(odds, 3)
        rec["p_fisher"] = p
        rec["p_bonferroni"] = min(1.0, p * N_MOTIFS)
        # Woolf 95% CI on log2 OR (Haldane-Anscombe 0.5 correction against zero cells)
        A, B, C, D = a + 0.5, b + 0.5, c_ + 0.5, d + 0.5
        lor = math.log((A * D) / (B * C))
        se = math.sqrt(1 / A + 1 / B + 1 / C + 1 / D)
        rec["log2_or"] = round(lor / math.log(2), 3)
        rec["log2_or_lo"] = round((lor - 1.96 * se) / math.log(2), 3)
        rec["log2_or_hi"] = round((lor + 1.96 * se) / math.log(2), 3)
        # enriched only if the whole CI sits above 1x (log2 OR > 0) after correction
        rec["enriched"] = bool(rec["log2_or_lo"] > 0 and rec["p_bonferroni"] < 0.05)
        rows.append(rec)
    spec = pd.DataFrame(rows).sort_values(["own_class", "motif"])
    spec.to_csv(HERE / "motif_specificity_matrix.csv", index=False)

    # ---- 2. per-protein motif status (incl. Mixed/Unassigned) ----
    g = (hits.drop_duplicates(["sequence_name", "motif"])
             .groupby("sequence_name")
             .agg(motifs=("motif", lambda s: ";".join(sorted(set(s)))),
                  n_motifs=("motif", "nunique")))
    prot = corpus.rename(columns={"Entry": "sequence_name"}).merge(
        g, on="sequence_name", how="left")
    prot["motifs"] = prot["motifs"].fillna("")
    prot["n_motifs"] = prot["n_motifs"].fillna(0).astype(int)
    prot["n_own"] = [sum(1 for m in ms.split(";") if m and m.split("-")[0] == c)
                     for ms, c in zip(prot["motifs"], prot["class"])]
    prot["n_foreign"] = prot["n_motifs"] - prot["n_own"]
    prot.rename(columns={"sequence_name": "Entry"}).to_csv(
        HERE / "fimo_cross_protein_motifs.csv", index=False)

    # ---- 3. report ----
    L = []
    L.append(f"FIMO cross-scan: 32 motifs x {len(corpus):,} proteins (p<1e-4)")
    L.append(f"  hits: {len(hits):,} over {hits.sequence_name.nunique():,} proteins\n")
    L.append("SPECIFICITY  (rate = % of that class's proteins carrying the motif)")
    L.append("  ratio_vs_bg  = own_pct / pooled-background_pct   <- primary")
    L.append("  ratio_vs_max = own_pct / worst single off-class   (blind to diffuse motifs)")
    L.append("  log2_or CI   = 95% Woolf CI; 'enriched' = CI entirely > 0 AND Bonferroni p < 0.05\n")
    show = ["motif", "consensus", "own_pct", "background_pct", "ratio_vs_bg",
            "ratio_vs_max", "worst_offclass", "log2_or", "log2_or_lo", "log2_or_hi",
            "p_bonferroni", "enriched"]
    L.append(spec[show].to_string(index=False))

    # the two ratios disagree where a motif is diffusely common rather than
    # concentrated in one rival class -- exactly the case ratio_vs_max misses
    dis = spec[(spec.ratio_vs_max >= 2) & (spec.ratio_vs_bg < 2)]
    if len(dis):
        L.append(f"\nDIFFUSE motifs -- pass vs-max but fail vs-background ({len(dis)}): "
                 "common across ALL other classes rather than one rival")
        L.append(dis[["motif", "consensus", "own_pct", "background_pct",
                      "ratio_vs_bg", "ratio_vs_max"]].to_string(index=False))

    weak = spec[spec.ratio_vs_bg < 2].sort_values("ratio_vs_bg")
    L.append(f"\nWEAK: <2x vs pooled background ({len(weak)} of {len(spec)}):")
    if len(weak):
        L.append(weak[["motif", "consensus", "own_pct", "background_pct",
                       "ratio_vs_bg", "log2_or_lo", "log2_or_hi", "enriched"]].to_string(index=False))
    depl = spec[spec.log2_or_hi < 0]
    L.append(f"\nDEPLETED in own class (CI entirely below 1x) -- misassigned to this class: {len(depl)}")
    if len(depl):
        L.append(depl[["motif", "consensus", "own_pct", "background_pct", "worst_offclass",
                       "log2_or", "log2_or_lo", "log2_or_hi"]].to_string(index=False))
    L.append("\nMIXED / UNASSIGNED (never scanned by the diagonal)")
    for c in ("Mixed", "Unassigned"):
        s = prot[prot["class"] == c]
        if not len(s):
            continue
        L.append(f"  {c}: {len(s):,} proteins | {(s.n_motifs > 0).sum():,} carry >=1 motif "
                 f"({100*(s.n_motifs > 0).mean():.1f}%)")
        top = (s[s.n_motifs > 0]["motifs"].str.split(";").explode()
                .value_counts().head(5))
        for m, n in top.items():
            L.append(f"      {m:12} {n:6,}")
    txt = "\n".join(L)
    (HERE / "fimo_cross_summary.txt").write_text(txt + "\n")
    print(txt)
    print("\n[wrote] motif_specificity_matrix.csv + fimo_cross_protein_motifs.csv + fimo_cross_summary.txt")


if __name__ == "__main__":
    main()
