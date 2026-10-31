#!/usr/bin/env python3
"""
SAM MTase UMAP by predicted substrate type -- ProstT5 version (datamapplot).

Same interactive datamapplot styling as sam_mtase_substrate_umap.py (clickable
legend toggles, search, hover, click-to-UniProt), but on the ProstT5 deployment
results instead of ProtT5:

  - UMAP coordinates + HDBSCAN clusters are read from the precomputed
    sam_mtase_umap_prostt5_coords.csv (global_reumap_prostt5.py output) -- no
    recompute.
  - Per-point substrate membership is reconstructed from the per-label scores
    and the frozen deployment thresholds (score >= threshold == a positive
    call), exactly reproducing the deployment call_* columns:
        one label positive    -> that substrate
        two or more positive  -> "Mixed"
        none positive         -> "None"

Outputs (sam_umap_output/):
  sam_mtase_umap_prostt5_static_datamap.png
  sam_mtase_umap_prostt5_interactive.html   (datamapplot; replaces the plotly one)

Env: UMAP (datamapplot + matplotlib).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import datamapplot
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).parent.resolve()
OUTPUT_DIR = SCRIPT_DIR
COORDS_CSV = OUTPUT_DIR / "sam_mtase_umap_prostt5_coords.csv"
THRESHOLDS_JSON = SCRIPT_DIR.parent / "2_training" / "deployment_model" / "thresholds.json"

LABELS = ["RNA", "DNA", "Protein", "Metabolite"]
SUBSTRATE_COLORS = {
    "RNA": "#cb181d", "DNA": "#6a51a3", "Protein": "#238b45",
    "Metabolite": "#e6550d", "Mixed": "#636363", "None": "#cccccc",
}


def get_predicted_substrate(membership):
    if membership == ["None"]:
        return "None"
    return membership[0] if len(membership) == 1 else "Mixed"


def create_hover_text(df):
    out = []
    for _, r in df.iterrows():
        out.append(
            f"Entry: {r['Entry']}\n"
            f"Organism: {r['Organism']}\n"
            f"Predicted: {' + '.join(r['Membership'])}\n"
            f"HDBSCAN cluster: {int(r['cluster'])}\n"
            f"---\n"
            f"RNA: {float(r['score_RNA']):.3f}\n"
            f"DNA: {float(r['score_DNA']):.3f}\n"
            f"Protein: {float(r['score_Protein']):.3f}\n"
            f"Metabolite: {float(r['score_Metabolite']):.3f}")
    return np.array(out)


def get_custom_legend(color_mapping):
    custom_css = """
    .row { display:flex; align-items:center; user-select:none; padding:2px 4px; border-radius:4px; }
    .legend-row { cursor:pointer; }
    .legend-row:hover { background:#00000010; }
    .legend-row.off { opacity:0.35; text-decoration:line-through; }
    .legend-key { cursor:default; }
    .legend-note { color:#888888; font-size:8.5pt; margin-left:4px; }
    .box { height:12px; width:12px; border-radius:3px; margin-right:8px; }
    #legend { position:absolute; bottom:0; right:0; margin:16px; padding:14px; border-radius:12px;
              z-index:2; background:#ffffffdd; font-family:Arial, sans-serif; font-size:11pt;
              box-shadow:2px 3px 10px #aaaaaa44; }
    #legend .legend-hint { margin-top:8px; font-size:8.5pt; color:#888888; }
    #title-container { max-width:75%; }
    """
    html = '<div id="legend"><b>Substrate Type</b><br><br>\n'
    for substrate, color in color_mapping.items():
        if substrate == "Mixed":
            html += (f'    <div class="row legend-key">'
                     f'<div id="{substrate}" class="box" style="background-color:{color};"></div>'
                     f'{substrate}<span class="legend-note">(shown under each label)</span></div>\n')
        else:
            html += (f'    <div class="row legend-row" data-substrate="{substrate}">'
                     f'<div id="{substrate}" class="box" style="background-color:{color};"></div>'
                     f'{substrate}</div>\n')
    html += '    <div class="legend-hint">Click a type to show / hide</div>\n</div>'
    return custom_css, html


def build_legend_toggle_js(df):
    toggle_categories = [c for c in SUBSTRATE_COLORS if c != "Mixed"]
    category_indices = {c: [] for c in toggle_categories}
    for i, members in enumerate(df["Membership"].tolist()):
        for label in members:
            category_indices.setdefault(label, []).append(i)
    js = """
(function () {
  const categoryIndices = __CATEGORY_INDICES__;
  const SELECTION_KIND = "substrate-legend";
  const visible = {};
  Object.keys(categoryIndices).forEach(c => { visible[c] = true; });
  function applyFilter() {
    if (Object.keys(visible).every(c => visible[c])) { datamap.removeSelection(SELECTION_KIND); return; }
    let union = [];
    for (const cat in categoryIndices) { if (visible[cat]) { union = union.concat(categoryIndices[cat]); } }
    if (union.length === 0) { union = [-1]; }
    datamap.addSelection(union, SELECTION_KIND);
  }
  function wireLegend() {
    document.querySelectorAll('#legend .legend-row').forEach(row => {
      row.addEventListener('click', () => {
        const cat = row.getAttribute('data-substrate');
        if (!(cat in visible)) { return; }
        visible[cat] = !visible[cat];
        row.classList.toggle('off', !visible[cat]);
        applyFilter();
      });
    });
  }
  function waitForReady() {
    if (typeof datamap !== 'undefined' && datamap.pointLayer &&
        document.querySelector('#legend .legend-row')) { wireLegend(); }
    else { setTimeout(waitForReady, 200); }
  }
  waitForReady();
})();
"""
    return js.replace("__CATEGORY_INDICES__", json.dumps(category_indices))


def create_interactive_plot(df, output_path):
    data_map = df[["UMAP_1", "UMAP_2"]].to_numpy().astype("float32")
    hover_data = create_hover_text(df)
    labels = np.where(df["cluster"].to_numpy() == -1, "Unlabelled",
                      np.array([f"Cluster {c}" for c in df["cluster"].to_numpy()]))
    color_array = df["Substrate"].map(SUBSTRATE_COLORS).fillna("#636363").values
    custom_css, custom_html = get_custom_legend(SUBSTRATE_COLORS)
    custom_js = build_legend_toggle_js(df)
    print("Creating interactive datamapplot...", flush=True)
    plot = datamapplot.create_interactive_plot(
        data_map, labels,
        font_family="Arial",
        title="SAM Methyltransferases by Predicted Substrate Type (ProstT5)",
        sub_title=f"UMAP projection of {len(df):,} tree-of-life SAM MTases (structure-aware embeddings)",
        hover_text=hover_data,
        enable_search=True,
        marker_color_array=color_array,
        custom_css=custom_css, custom_html=custom_html, custom_js=custom_js,
        marker_size_array=np.ones(len(df)) * 5,
        extra_point_data=df[["Entry"]].copy(),
        on_click="window.open(`https://www.uniprot.org/uniprotkb/{Entry}/entry`)",
    )
    plot.save(str(output_path))
    print(f"Saved interactive plot to: {output_path}", flush=True)


def create_static_plot(df, output_path):
    fig, ax = plt.subplots(figsize=(14, 10))
    draw_order = ["None", "Mixed", "Metabolite", "Protein", "DNA", "RNA"]
    dim = {"None": (0.25, 5), "Mixed": (0.40, 6)}
    handles = {}
    for k, substrate in enumerate(draw_order):
        mask = df["Substrate"] == substrate
        if mask.sum() == 0:
            continue
        alpha, size = dim.get(substrate, (0.7, 8))
        handles[substrate] = ax.scatter(
            df.loc[mask, "UMAP_1"], df.loc[mask, "UMAP_2"], c=SUBSTRATE_COLORS[substrate],
            label=f"{substrate} (n={mask.sum():,})", alpha=alpha, s=size,
            edgecolors="none", zorder=2 + k)
    ax.set_xlabel("UMAP 1", fontsize=12); ax.set_ylabel("UMAP 2", fontsize=12)
    legend_order = [s for s in SUBSTRATE_COLORS if s in handles]
    ax.legend([handles[s] for s in legend_order], [handles[s].get_label() for s in legend_order],
              loc="upper right", fontsize=10, framealpha=0.9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved static plot to: {output_path}", flush=True)


def main():
    df = pd.read_csv(COORDS_CSV)
    thr = json.load(open(THRESHOLDS_JSON))["thresholds"]
    print(f"Loaded {len(df):,} points; thresholds = {thr}", flush=True)

    # Reconstruct per-point membership from scores + frozen thresholds.
    def membership(row):
        called = [l for l in LABELS if float(row[f"score_{l}"]) >= thr[l]]
        return called if called else ["None"]
    df["Membership"] = df.apply(membership, axis=1)
    df["Substrate"] = df["Membership"].apply(get_predicted_substrate)
    print("Substrate distribution:", df["Substrate"].value_counts().to_dict(), flush=True)

    create_static_plot(df, OUTPUT_DIR / "sam_mtase_umap_prostt5_static_datamap.png")
    create_interactive_plot(df, OUTPUT_DIR / "sam_mtase_umap_prostt5_interactive.html")
    print("done", flush=True)


if __name__ == "__main__":
    main()
