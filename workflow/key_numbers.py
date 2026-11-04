#!/usr/bin/env python3
"""The numbers of a release that the Results quote, in one file.

  python workflow/key_numbers.py [--release <dir>]       (default: $RELEASE_DIR, else the newest release in $UNIPROT_PULL_DIR)

Writes <release>/key_numbers.md, read from the outputs of the other stages (never recomputed from raw data, so every number can be traced to a file, named in the heading of its section).
A section whose input is missing is listed as "not available" with the reason: the stage that makes it did not run or failed. Regenerated at the end of the refresh chain, so after a
new UniProt release this file is the checklist for updating the draft.
"""
import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "1_clustering"))
SECTIONS = []


def section(title):
    def deco(fn):
        SECTIONS.append((title, fn))
        return fn
    return deco


def pct(x, n=1):
    return f"{100 * x:.{n}f}%"


def fmt_p(p):
    """'< 1e-4' / '= 3.4e-03' / '= 0.050': the sign included, so a caller writes `p {fmt_p(p)}`."""
    return "< 1e-4" if p < 1e-4 else f"= {p:.1e}" if p < 1e-2 else f"= {p:.3f}"


@section("Data  (pull_manifest.json, labels_SAM.csv, proteome_summary.csv)")
def data(rel):
    m = json.load(open(rel / "pull_manifest.json"))
    d = m["datasets"]
    out = [f"UniProt release {m['uniprot_release']} ({m.get('uniprot_release_date', '')})",
           f"reviewed InterPro IPR029063 entries: {d['sam_mtases']['rows']:,} ({d['sam_mtases']['rows_le_length_view']:,} up to 2,000 aa)",
           f"reviewed enzymes with a catalytic-activity annotation (UP): {d['up_enzymes']['rows']:,} ({d['up_enzymes']['rows_le_length_view']:,} up to 2,000 aa)"]
    lab = pd.read_csv(rel / "labels_SAM.csv", low_memory=False)
    if "in_benchmark" in lab:
        b = lab[lab.in_benchmark.astype(bool)]
        out.append(f"benchmark set (cluster-disjoint cross-validation): {len(b):,} proteins; labels " + ", ".join(f"{c} {int(b[c].sum()):,}" for c in ["Metabolite", "Protein", "RNA", "DNA"] if c in b))
    ps = pd.read_csv(rel / "proteome_summary.csv")
    ok = ps[ps.status == "ok"]
    out.append(f"reference proteomes: {len(ok):,}; proteins in them: {int(ok.n_all.sum()):,}; with IPR029063: {int(ok.n_ipr_all.sum()):,} ({int(ok.n_ipr_le.sum()):,} up to 2,000 aa)")
    return out


@section("Benchmark  (benchmark/sam_classifier_comparison_summary.csv, stats/)")
def bench(rel):
    s = pd.read_csv(rel / "benchmark" / "sam_classifier_comparison_summary.csv")
    s = s[(s.split_mode == "cluster") & (s.dataset == "full")].sort_values("mcc_mean", ascending=False)
    cfg = json.load(open(REPO / "refresh_config.json"))
    ctrl = cfg.get("control", "MLP_ProtT5")
    out = [f"deployed (control) model: {ctrl}", "cluster-disjoint MCC, mean ± SD over 16 iterations: " + "; ".join(f"{r.classifier} {r.mcc_mean:.3f} ± {r.mcc_std:.3f}" for r in s.itertuples())]
    f = rel / "benchmark" / "stats" / "friedman_fig2A.txt"
    if f.exists():
        m = re.search(r"chi-?square\((\d+)\)\s*=\s*([\d.]+),\s*p\s*=\s*([\d.e+-]+)", f.read_text()) or re.search(r"chi2\((\d+)\)=([\d.]+).*?p=([\d.e+-]+)", f.read_text())
        if m:
            out.append(f"Friedman test, Fig 2A methods: chi-square({m.group(1)}) = {m.group(2)}, p = {m.group(3)}")
    t = rel / "benchmark" / "stats" / "table_S4_wilcoxon_holm.csv"
    if t.exists():
        t4 = pd.read_csv(t)
        sig = t4[t4["holm_p"] < 0.05]
        better, worse = sig[sig.mean_diff_vs_control > 0], sig[sig.mean_diff_vs_control < 0]
        out.append(f"Table S4 (each method vs {ctrl}, Wilcoxon, Holm over {len(t4)} comparisons): significantly worse: " + (", ".join(worse.method) or "none")
                   + "; significantly BETTER: " + (", ".join(f"{r.method} (MCC {r.mean_mcc:.3f}, wins {r.iterations_better}/{r.n_iterations}, Holm p {fmt_p(r.holm_p)})" for r in better.itertuples()) or "none"))
    return out


@section("MTase-trained vs UP-trained  (sam_vs_up/sam_vs_up_paired_tests.csv)")
def sam_up(rel):
    t = pd.read_csv(rel / "sam_vs_up" / "sam_vs_up_paired_tests.csv")
    t = t[(t.comparison == "ProtT5_UP minus ProtT5_SAM") & (t.metric == "mcc")]
    out = []
    for r in t.itertuples():
        out.append(f"{r.label}: MTase-trained {r.mean_b:.3f}, UP-trained {r.mean_a:.3f}, difference {r.mean_diff:+.3f} (UP wins {r.splits_a_wins}); Wilcoxon p {fmt_p(r.p_wilcoxon)}, Nadeau-Bengio p {fmt_p(r.p_nadeau_bengio)}")
    return out


@section("Proteome predictions  (predictions_all_proteomes.csv.gz)")
def preds(rel):
    p = pd.read_csv(rel / "predictions_all_proteomes.csv.gz", low_memory=False)
    calls = p[[f"call_{c}" for c in ["Metabolite", "Protein", "RNA", "DNA"]]].astype(int)
    n = len(p)
    k = calls.sum(axis=1)
    out = [f"{n:,} proteins from {p.proteome.nunique():,} proteomes (up to 2,000 aa)",
           "called (score ≥ 0.5): " + ", ".join(f"{c} {pct(calls['call_' + c].mean())}" for c in ["Metabolite", "Protein", "RNA", "DNA"]),
           f"exactly one call {pct((k == 1).mean())}; more than one {pct((k > 1).mean())}; none {pct((k == 0).mean())}"]
    if "in_domain" in p:
        out.append(f"inside the applicability domain: {pct(p.in_domain.astype(int).mean())}")
    return out


@section("Phylogenetic trends  (fig4/fig4_statistics.csv, fig5/fig5_scaling_numbers.csv, fig5/lineage_statistics.txt)")
def trends(rel):
    s = pd.read_csv(rel / "fig4" / "fig4_statistics.csv")
    out = []
    for r in s.itertuples():
        if pd.isna(r.n_b):      # an omnibus (Kruskal-Wallis) row
            out.append(f"Fig 4{r.panel} {r.measure}: {r.group_a} (n={int(r.n_a)}); {r.test} p {fmt_p(r.p)}")
        else:
            out.append(f"Fig 4{r.panel} {r.measure}: {r.group_a} (n={int(r.n_a)}, median {r.median_a:.4g}) vs {r.group_b} (n={int(r.n_b)}, median {r.median_b:.4g}); {r.test} p {fmt_p(r.p)}")
    nf = rel / "fig5" / "fig5_scaling_numbers.csv"
    if nf.exists():      # Fig 5 as power-law fits n = lambda * g^alpha (figures/Fig5_proteome_scaling/scaling_loglog.py)
        n = pd.read_csv(nf)
        P = n[n.analysis == "primary"]
        for pl, g in P.groupby("panel", sort=True):
            q = lambda grp, qty: g[(g.group == grp) & (g.quantity == qty)].iloc[0]
            lin, rest = [x for x in g.group.unique() if x not in ("lineage - rest", "lineage / rest")][:2]
            a1, a2, d, f = q(lin, "alpha"), q(rest, "alpha"), q("lineage - rest", "delta_alpha"), q("lineage / rest", "fold_at_median_lineage_size")
            lg = n[(n.panel == pl) & (n.analysis == "legacy_ols_slope_ratio")].estimate.iloc[0]
            out.append(f"Fig 5{pl} {g.measure.iloc[0]}: alpha {lin} {a1.estimate:.2f} [{a1.ci_lo:.2f}, {a1.ci_hi:.2f}], {rest} {a2.estimate:.2f} [{a2.ci_lo:.2f}, {a2.ci_hi:.2f}]; "
                       f"delta alpha {d.estimate:+.2f} [{d.ci_lo:+.2f}, {d.ci_hi:+.2f}]; level at equal size x{f.estimate:.2f} [{f.ci_lo:.2f}, {f.ci_hi:.2f}] ({f.note.split(';')[0]}); "
                       f"(old OLS slope ratio {lg:.2f}: not a rate, do not quote)")
    lt = rel / "fig5" / "lineage_statistics.txt"
    if lt.exists():
        out += lt.read_text().strip().splitlines()
    return out


@section("Map and clusters  (umap/umap_report.txt, umap/cluster_selection_report.txt)")
def umap(rel):
    out = [l for l in (rel / "umap" / "umap_report.txt").read_text().splitlines() if l.strip()][:5]
    out += [l for l in (rel / "umap" / "cluster_selection_report.txt").read_text().splitlines()[:6] if l.strip()]
    return out


@section("Dark clusters  (dark_annotation/table_dark_clusters.csv, horizyn_run/horizyn_predict_full.log, fig8/fig8_clusters.csv)")
def dark(rel):
    t = pd.read_csv(rel / "dark_annotation" / "table_dark_clusters.csv")
    t["tier"] = t.clean_confidence.str.split().str[0]
    out = [f"{len(t)} dark clusters: " + ", ".join(f"{k} {v}" for k, v in t.lineage.value_counts().items()) + f"; {int(t.n.sum()):,} proteins"]
    for lin, g in t.groupby("lineage"):
        out.append(f"{lin}: CLEAN reliable {int((g.tier == 'reliable').sum())}, marginal {int((g.tier == 'marginal').sum())}, unreliable {int((g.tier == 'unreliable').sum())}; called a methyltransferase (EC 2.1.1.*) {int((g.clean_mtase == 'yes').sum())} of {len(g)}")
    log = rel / "dark_annotation" / "horizyn_run" / "horizyn_predict_full.log"
    if log.exists():
        txt = log.read_text()
        m = re.search(r"median best-SAM-methylation rank: (\d+) of (\d+)", txt)
        m2 = re.search(r"proteins with >=1 SAM-methylation reaction in their own top-100: (\d+) / (\d+)", txt)
        if m and m2:
            out.append(f"Horizyn: median rank of the best SAM-methylation reaction {m.group(1)} of {int(m.group(2)):,}; {m2.group(1)} of {m2.group(2)} proteins have one in their top 100")
    f8 = rel / "fig8" / "fig8_clusters.csv"
    if f8.exists():
        out.append("Fig 8 clusters (CLEAN reliable + methyltransferase, EC inside Horizyn's dominant family): " + ", ".join(f"c{c}" for c in pd.read_csv(f8).cluster))
    return out


@section("Motifs  (motif_analysis/motif_specificity_matrix.csv, motif_status_all.csv, structure_search/motif_placements.csv)")
def motifs(rel):
    d = rel / "motif_analysis"
    s = pd.read_csv(d / "motif_specificity_matrix.csv")
    kept = s[s.ratio_vs_bg >= 2]
    out = [f"{len(kept)} of {len(s)} motifs are at least 2-fold enriched in their class; by class: " + ", ".join(f"{c} {int((kept.own_class == c).sum())}" for c in ["DNA", "RNA", "Protein", "Metabolite"]),
           "below 2-fold: " + ", ".join(f"{r.motif} {r.consensus} ({r.ratio_vs_bg:.1f}x)" for r in s[s.ratio_vs_bg < 2].itertuples())]
    ms = pd.read_csv(d / "motif_status_all.csv")
    out.append(f"proteins without any motif of their class: {int((~ms.has_motif).sum()):,} of {len(ms):,} ({pct((~ms.has_motif).mean())}); " + ", ".join(f"{c} {pct((~g.has_motif).mean())}" for c, g in ms.groupby("class")))
    pl = pd.read_csv(d / "structure_search" / "motif_placements.csv")
    out.append(f"structure of every motif: {len(pl)} of {len(kept)} found; placement " + ", ".join(f"{k} {v}" for k, v in pl.placement.value_counts().items()))
    weak = pl[(~pl.in_class) | pl.substrate_kind_fallback]
    if len(weak):
        out.append("weak structure picks (other class or other substrate kind): " + ", ".join(weak.motif))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release")
    a = ap.parse_args()
    from mmseqs_clusters import newest_release
    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    try:
        tag = "UniProt release " + json.load(open(rel / "pull_manifest.json"))["uniprot_release"]
    except Exception:
        tag = rel.name
    md = [f"# Key numbers of {tag}", "", "Read from the outputs of the stages (the file of each section is named in its heading); regenerate with `python workflow/key_numbers.py`.", ""]
    for title, fn in SECTIONS:
        md.append(f"## {title}")
        try:
            md += [f"- {l}" for l in fn(rel)]
        except Exception as e:      # a missing input: say which
            md.append(f"- not available: {type(e).__name__}: {e}")
        md.append("")
    (rel / "key_numbers.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
