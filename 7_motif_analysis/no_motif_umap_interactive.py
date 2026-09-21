#!/usr/bin/env python3
"""
Interactive motif-dark UMAP (datamapplot) with an A/B view toggle.

Same datamapplot infrastructure as UMAP/sam_mtase_substrate_umap_prostt5.py
(hover, search, click-to-UniProt), but the point set is highlighted by STREME
motif coverage instead of a per-substrate legend toggle:

  Every point is colored by substrate. A "view" then SELECTS a subset (which
  stays in full color) and dims everything else to grey -- reusing datamapplot's
  addSelection / removeSelection, exactly the mechanism the substrate script uses.

  View A  "Motif-dark families": points in UMAP (HDBSCAN) clusters whose every
          labeled member is no-motif (11 clusters).
  View B  "All no-motif":       every STREME-negative protein (11,112), which are
          mostly scattered inside motif families.

Two-level labeling (see no_motif_umap.py): mmseqs 50%-id rep motif status
propagated to all members, then UMAP-cluster purity.

Env: UMAP (datamapplot). Output: no_motif_umap_interactive.html
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import datamapplot

MA = Path(__file__).resolve().parent
COORDS = Path("/home/nroethler/Code/journal_article/repo/5_umap/sam_mtase_umap_prostt5_coords.csv")
CLUST = MA / "clust_cluster.tsv"
NOMOTIF = MA / "attribution/no_motif_manifest.csv"
OUT_HTML = MA / "no_motif_umap_interactive.html"

SUBSTRATE_COLORS = {"RNA": "#cb181d", "DNA": "#6a51a3", "Protein": "#238b45",
                    "Other": "#e6550d", "Mixed": "#636363", "Unassigned": "#dddddd"}


def load_labeled():
    umap = pd.read_csv(COORDS)
    clu = pd.read_csv(CLUST, sep="\t", header=None, names=["rep", "member"])
    nm = pd.read_csv(NOMOTIF)
    no_motif_reps = set(nm["Entry"])
    clu["rep_bare"] = clu["rep"].str.split("|").str[0]
    clu["member_bare"] = clu["member"].str.split("|").str[0]
    clu["status"] = np.where(clu["rep_bare"].isin(no_motif_reps), "no-motif", "motif")
    umap["status"] = umap["Entry"].map(dict(zip(clu["member_bare"], clu["status"])))
    return umap.reset_index(drop=True)


def pure_clusters(umap):
    keep = umap[umap["cluster"] != -1]
    g = keep.groupby("cluster")["status"]
    n_motif = g.apply(lambda s: (s == "motif").sum())
    n_no = g.apply(lambda s: (s == "no-motif").sum())
    return set(n_motif[(n_motif == 0) & (n_no > 0)].index)


def hover_text(df):
    out = []
    for _, r in df.iterrows():
        st = r["status"] if isinstance(r["status"], str) else "not in motif corpus"
        out.append(f"Entry: {r['Entry']}\nOrganism: {r['Organism']}\n"
                   f"Substrate: {r['Substrate']}\nMotif status: {st}\n"
                   f"HDBSCAN cluster: {int(r['cluster'])}")
    return np.array(out)


CUSTOM_CSS = """
#viewbox { position:absolute; top:0; left:0; margin:16px; padding:14px; border-radius:12px;
           z-index:2; background:#ffffffdd; font-family:Arial, sans-serif; font-size:11pt;
           box-shadow:2px 3px 10px #aaaaaa44; max-width:260px; }
#viewbox b { font-size:12pt; }
.viewbtn { display:block; width:100%; text-align:left; margin:6px 0 0; padding:8px 10px;
           border:1px solid #ccc; border-radius:8px; background:#f7f7f7; cursor:pointer;
           font-family:Arial; font-size:10.5pt; }
.viewbtn:hover { background:#eee; }
.viewbtn.active { background:#2c7fb8; color:#fff; border-color:#2c7fb8; font-weight:bold; }
.viewnote { color:#888; font-size:8.5pt; margin-top:4px; }
#leg { position:absolute; bottom:0; right:0; margin:16px; padding:12px; border-radius:12px;
       z-index:2; background:#ffffffdd; font-family:Arial; font-size:10.5pt;
       box-shadow:2px 3px 10px #aaaaaa44; }
#leg .row { display:flex; align-items:center; padding:2px 4px; }
#leg .box { height:12px; width:12px; border-radius:3px; margin-right:8px; }
#title-container { max-width:70%; }
"""


def custom_html():
    h = ('<div id="viewbox"><b>Motif coverage</b>'
         '<div class="viewnote">Colored = shown; grey = dimmed.</div>'
         '<button class="viewbtn" id="btnA" data-view="A">A &nbsp;Motif-dark families'
         '<div class="viewnote">HDBSCAN clusters with no motif in any member</div></button>'
         '<button class="viewbtn active" id="btnB" data-view="B">B &nbsp;All no-motif proteins'
         '<div class="viewnote">every STREME-negative protein</div></button>'
         '<button class="viewbtn" id="btnAll" data-view="ALL">Show all</button></div>')
    h += '<div id="leg"><b>Substrate</b><br>'
    for s in ["RNA", "DNA", "Protein", "Other", "Mixed", "Unassigned"]:
        h += (f'<div class="row"><div class="box" style="background:{SUBSTRATE_COLORS[s]}"></div>{s}</div>')
    h += '</div>'
    return h


def custom_js(idx_A, idx_B):
    js = """
(function () {
  const IDX = {"A": __IDXA__, "B": __IDXB__};
  const KIND = "motif-view";
  function applyView(v) {
    document.querySelectorAll('.viewbtn').forEach(b =>
      b.classList.toggle('active', b.getAttribute('data-view') === v));
    if (v === 'ALL') { datamap.removeSelection(KIND); return; }
    datamap.removeSelection(KIND);
    datamap.addSelection(IDX[v], KIND);
  }
  function wire() {
    document.querySelectorAll('.viewbtn').forEach(b =>
      b.addEventListener('click', () => applyView(b.getAttribute('data-view'))));
    applyView('B');   // default landing view
  }
  function ready() {
    if (typeof datamap !== 'undefined' && datamap.pointLayer &&
        document.querySelector('.viewbtn')) { wire(); }
    else { setTimeout(ready, 200); }
  }
  ready();
})();
"""
    return js.replace("__IDXA__", json.dumps(idx_A)).replace("__IDXB__", json.dumps(idx_B))


def main():
    df = load_labeled()
    pure = pure_clusters(df)
    idx_A = np.where(df["cluster"].isin(pure).to_numpy())[0].tolist()
    idx_B = np.where((df["status"] == "no-motif").to_numpy())[0].tolist()
    print(f"points {len(df):,} | view A {len(idx_A):,} (in {len(pure)} pure clusters) | "
          f"view B {len(idx_B):,} no-motif", flush=True)

    data_map = df[["UMAP_1", "UMAP_2"]].to_numpy().astype("float32")
    labels = np.where(df["cluster"].to_numpy() == -1, "Unlabelled",
                      np.array([f"Cluster {c}" for c in df["cluster"].to_numpy()]))
    color_array = df["Substrate"].map(SUBSTRATE_COLORS).fillna("#bbbbbb").values

    plot = datamapplot.create_interactive_plot(
        data_map, labels,
        font_family="Arial",
        title="SAM-MTase landscape — STREME motif coverage",
        sub_title="Toggle A (entirely motif-dark families) / B (all no-motif proteins). Colored by substrate.",
        hover_text=hover_text(df),
        enable_search=True,
        marker_color_array=color_array,
        custom_css=CUSTOM_CSS, custom_html=custom_html(), custom_js=custom_js(idx_A, idx_B),
        marker_size_array=np.ones(len(df)) * 5,
        extra_point_data=df[["Entry"]].copy(),
        on_click="window.open(`https://www.uniprot.org/uniprotkb/{Entry}/entry`)",
    )
    plot.save(str(OUT_HTML))
    print(f"[wrote] {OUT_HTML}", flush=True)


if __name__ == "__main__":
    main()
