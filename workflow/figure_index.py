#!/usr/bin/env python3
"""An index of every figure and table of a release, for the review before anything is pushed.

  python workflow/figure_index.py [--release <dir>]       (default: $RELEASE_DIR, else the newest release in $UNIPROT_PULL_DIR)

Writes <release>/FIGURES.md (a table: item, file, stage that makes it, when it was made, what to look at) and <release>/FIGURES.html (the same as a gallery of the images, to page through; one self-contained file with downscaled copies of the images, so it can be opened or sent anywhere).
A file that is missing is listed as MISSING: the stage that makes it did not run or failed.
"""
import argparse
import base64
import datetime
import html
import io
import os
import shutil
import sys
from pathlib import Path

from PIL import Image

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "1_clustering"))

# (label, file relative to the release, stage of workflow/run_refresh.sh, what to look at)
ITEMS = [
    ("Fig 2  A: representations and search tools; B: MTase-trained vs UP-trained", "benchmark/Fig2.png", "fig2", "A: cluster-split MCC, stars = Wilcoxon-Holm vs MLP_ProtT5, Friedman box; B: Wilcoxon-Holm brackets; statistics sentences in benchmark/Fig2_caption.txt"),
    ("Table 1  metrics", "benchmark/metrics_table.png", "report", "numbers quoted in the Results"),
    ("Fig S1  train/validation loss, benchmark MLPs", "benchmark/train_val_loss_by_embedding.png", "report", "supplementary figures are numbered in the order the draft first cites them: this is cited first (Methods, loss settings)"),
    ("Fig S2  train/validation loss, UP-trained MLP", "sam_vs_up/up_train_val_loss.png", "sam_vs_up_report", "cited after Fig S1"),
    ("(not numbered)  per-label loss, benchmark MLPs", "benchmark/train_val_loss_by_label.png", "report", "per-label version of Fig S1; not cited in the draft"),
    ("Fig 3  tree of life", "tree/fig3_tree.png", "tree", "rings, bars, legend; also .svg / .pdf"),
    ("Fig 4  trends (composite)", "fig4/Fig4_trends_combined.png", "fig4", "A: Protein share by cellularity (the difference vanished on the ProtT5 calls); B-D: violins; statistics in fig4/fig4_statistics.csv"),
    ("Fig 5  proteome size vs SAM-MTases, by clade", "fig5/fig5_by_clade.png", "fig5", "one power-law fit per clade (A, B Actinomycetota and other bacteria; C, D plants, fungi, metazoa, other eukaryotes), clade-bootstrap intervals (fig5/fig5_by_clade.csv, fig5/fig5_clade_model_comparison.csv)"),
    ("(not numbered)  Fig 5 previous version, plants vs other eukaryotes and Actinomycetota vs other bacteria", "fig5/fig5_scaling_loglog.png", "fig5", "the lineage-vs-rest fits (fig5/fig5_scaling_numbers.csv); the eukaryote rest group mixes clades at different levels"),
    ("Fig 6  global UMAP", "umap/sam_mtase_umap_prott5_static_datamap.png", "umap", "also the interactive umap/sam_mtase_umap_prott5_interactive.html"),
    ("Fig 7  cluster grids (composite)", "umap/fig7/Fig7_clusters_combined.png", "annotate", "A metazoan-protein, B vertebrate-DNA, C plant dark, D actinomycete dark clusters; Fig 8's clusters outlined"),
    ("Fig 8  Horizyn reaction profiles", "fig8/Fig8_reaction_profiles.png", "annotate", "clusters chosen by rule; check the panel titles and RHEA numbers"),
    ("Fig 9  motifs at known sites", "fig9_10/Fig9_Known.png", "fig9_10", "closest SAM + substrate motif per class; check the zoom and the captions"),
    ("Fig 10  motifs at distal sites", "fig9_10/Fig10_Distal.png", "fig9_10", "furthest distal motif per class"),
    ("Fig 11  no-motif UMAP", "motif_analysis/no_motif_umap.png", "no_motif_umap", "regions of the map without a motif"),
    ("Table 2  motifs", "motif_analysis/motif_table.png", "motif_structures", "24 motifs of at least 2-fold; Protein-2, 3, 6, 8 rest on weak structure picks"),
    ("Table S1  Actinomycetota dark clusters", "dark_annotation/tables/table_s1.png", "annotate", "biological-context notes and family descriptions still to be written"),
    ("Table S2  plant dark clusters", "dark_annotation/tables/table_s2.png", "annotate", "same"),
    ("Table S3  InterPro domain enrichment", "motif_analysis/domain_enrichment_table.png", "domain_table", ""),
    ("Numbers  lineage comparisons", "fig5/lineage_statistics.txt", "fig5", "the lineage statistics quoted in the Results"),
    ("Numbers  Fig 5 fits", "fig5/fig5_scaling_numbers.csv", "fig5", ""),
    ("Numbers  Fig 5 clade fits", "fig5/fig5_by_clade.csv", "fig5", ""),
    ("Numbers  Fig 4 statistics", "fig4/fig4_statistics.csv", "fig4", ""),
]
# the files for the article draft: full-size copies under these names in <release>/article_figures/ (the same files as listed in the index)
EXPORT = {"benchmark/Fig2.png": "Fig2.png", "benchmark/Fig2_caption.txt": "Fig2_caption.txt", "benchmark/metrics_table.png": "Table1.png",
          "benchmark/train_val_loss_by_embedding.png": "FigS1.png", "sam_vs_up/up_train_val_loss.png": "FigS2.png",
          "benchmark/train_val_loss_by_label.png": "unnumbered_train_val_loss_by_label.png", "tree/fig3_tree.png": "Fig3.png", "tree/fig3_tree.svg": "Fig3.svg", "tree/fig3_tree.pdf": "Fig3.pdf",
          "fig4/Fig4_trends_combined.png": "Fig4.png", "fig5/fig5_by_clade.png": "Fig5.png", "fig5/fig5_scaling_loglog.png": "unnumbered_fig5_lineage_vs_rest.png", "umap/sam_mtase_umap_prott5_static_datamap.png": "Fig6.png",
          "umap/fig7/Fig7_clusters_combined.png": "Fig7.png", "fig8/Fig8_reaction_profiles.png": "Fig8.png", "fig9_10/Fig9_Known.png": "Fig9.png", "fig9_10/Fig10_Distal.png": "Fig10.png",
          "motif_analysis/no_motif_umap.png": "Fig11.png", "motif_analysis/motif_table.png": "Table2.png", "dark_annotation/tables/table_s1.png": "TableS1.png",
          "dark_annotation/tables/table_s2.png": "TableS2.png", "motif_analysis/domain_enrichment_table.png": "TableS3.png"}
IMG = (".png", ".jpg", ".jpeg")
MAX_W = 1500        # width of the embedded copies; the full-size files are the ones listed


def data_uri(path):
    """A downscaled copy of the image as a data: URI, so FIGURES.html shows its images wherever the file is opened (relative links break as soon as it is moved or opened remotely)."""
    im = Image.open(path).convert("RGB")
    if im.width > MAX_W:
        im = im.resize((MAX_W, round(im.height * MAX_W / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "JPEG", quality=88, optimize=True)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()



def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", default=None)
    a = ap.parse_args()
    from mmseqs_clusters import newest_release
    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    rows, miss = [], 0
    for label, f, stage, note in ITEMS:
        p = rel / f
        if p.exists():
            when = datetime.datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        else:
            when, miss = "MISSING", miss + 1
        rows.append((label, f, stage, when, note))
    exp = rel / "article_figures"
    exp.mkdir(exist_ok=True)
    for f, name in EXPORT.items():
        if (rel / f).exists():
            shutil.copyfile(rel / f, exp / name)
    md = [f"# Figures and tables of {rel.name}", "", f"Made by `workflow/figure_index.py`; paths are relative to `{rel}`. {miss} missing.", "",
          "| Item | File | Stage | Made | Look at |", "|---|---|---|---|---|"]
    md.insert(3, f"For the article draft, the full-size files are copied under their figure numbers into `{exp}/` (Fig2A.png, Fig2B.png, Table1.png, Fig3.png/.svg/.pdf, ...).\n")
    md += [f"| {l} | `{f}` | {s} | {w} | {n} |" for l, f, s, w, n in rows]
    (rel / "FIGURES.md").write_text("\n".join(md) + "\n")
    cards = []
    for l, f, s, w, n in rows:
        body = (f'<img src="{data_uri(rel / f)}" alt="{html.escape(l)}">' if f.lower().endswith(IMG) and w != "MISSING"
                else "")
        cards.append(f'<section><h3>{html.escape(l)}</h3><p class="m">{html.escape(f)} &middot; stage {html.escape(s)} &middot; {html.escape(w)}</p><p>{html.escape(n)}</p>{body}</section>')
    (rel / "FIGURES.html").write_text("<!doctype html><meta charset=utf-8><title>Figures</title><style>body{font:15px sans-serif;max-width:1100px;margin:2em auto;padding:0 1em}"
                                      "section{border-top:1px solid #ccc;padding:1em 0}img{max-width:100%}.m{color:#666;font-size:12px}</style><h1>Figures and tables of "
                                      f"{html.escape(rel.name)}</h1>" + "".join(cards))
    print("\n".join(md))
    print(f"\nwrote {rel/'FIGURES.md'} and {rel/'FIGURES.html'}; full-size copies for the draft in {exp}/")


if __name__ == "__main__":
    main()
