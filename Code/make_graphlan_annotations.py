#!/usr/bin/env python3
import os
import re
import argparse
from pathlib import Path
import pandas as pd

# ============ CONFIG (defaults; override via CLI) ============
PREDICTION_DIR = "prediction_files"
TREE_NEWICK = "taxonomy_tree_sanitized.nwk"            # your .nwk from the lineage step
ANNOT_OUT = "graphlan_annotations.txt"      # annotation file to feed graphlan_annotate.py

MODE = "percent"        # "percent" or "mean"
THRESH = 0.5            # threshold for "percent" mode
SCALE_RING_HEIGHT = False  # if True, ring height varies with value; else constant height
RING_LEVEL = 1
TITLE = "RNA-binding across species (SAM MTases)"
RING_LABEL_PERCENT = "RNA-binding % (≥{:.2f})".format(THRESH)
RING_LABEL_MEAN = "Mean RNA score"

# -------- utility: match the Newick escaping we used when writing the tree --------
def escape_newick_label(label: str) -> str:
    if label is None:
        return ""
    label = label.replace(" ", "_")
    label = re.sub(r"[,:;()]", "_", label)
    return label

# -------- simple linear gradient (light gray → blue) without matplotlib --------
def lerp(a, b, t):
    return int(round(a + (b - a) * t))

def hex_to_rgb(hexcolor: str):
    hexcolor = hexcolor.lstrip("#")
    return tuple(int(hexcolor[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_hex(rgb):
    return "#{:02X}{:02X}{:02X}".format(*rgb)

# Nice readable gradient: light gray → deep blue
C0 = "#e9eef3"  # low
C1 = "#084594"  # high

def value_to_hex(v: float, vmin=0.0, vmax=1.0) -> str:
    if v is None:
        return "#BDBDBD"  # neutral gray for missing
    # clamp
    if v < vmin: v = vmin
    if v > vmax: v = vmax
    t = 0.0 if vmax == vmin else (v - vmin) / (vmax - vmin)
    r0, g0, b0 = hex_to_rgb(C0)
    r1, g1, b1 = hex_to_rgb(C1)
    rgb = (lerp(r0, r1, t), lerp(g0, g1, t), lerp(b0, b1, t))
    return rgb_to_hex(rgb)

# -------- collect per-species values from prediction CSVs --------
def collect_species_values(pred_dir: str, mode: str, thresh: float):
    """
    Returns:
        species_val: dict {species_label_for_tree: value in [0,1]}
        species_n:   dict {species_label_for_tree: count of proteins aggregated}
    The species label is escaped to match the .nwk leaf labels.
    """
    species_records = {}  # species -> list of RNA floats
    files = 0
    for fname in os.listdir(pred_dir):
        if not fname.endswith(".csv"):
            continue
        fpath = os.path.join(pred_dir, fname)
        try:
            df = pd.read_csv(fpath, dtype=str, encoding_errors="ignore", on_bad_lines="skip")
        except Exception as e:
            print(f"⚠️ Could not read {fname}: {e}")
            continue
        files += 1
        # find columns
        org_col = next((c for c in ["Organism","organism","Species","species"] if c in df.columns), None)
        rna_col = next((c for c in ["RNA","rna","Rna"] if c in df.columns), None)
        if org_col is None or rna_col is None:
            print(f"⚠️ Missing columns in {fname} (found: {list(df.columns)[:6]} …)")
            continue

        # cast RNA to float; coerce errors to NaN
        rna = pd.to_numeric(df[rna_col], errors="coerce")
        org = df[org_col].astype(str)

        sub = pd.DataFrame({"org": org, "rna": rna}).dropna()
        if sub.empty: 
            continue

        # group per species
        for sp, grp in sub.groupby("org"):
            vals = grp["rna"].astype(float).tolist()
            if not vals:
                continue
            species_records.setdefault(sp, []).extend(vals)

    print(f"📄 Scanned {files} CSV file(s) under {pred_dir}")
    # aggregate
    species_val = {}
    species_n = {}
    for sp, vals in species_records.items():
        if not vals:
            continue
        if mode == "mean":
            v = float(pd.Series(vals).mean())
        else:  # percent
            v = float((pd.Series(vals) >= thresh).mean())
        species_val[escape_newick_label(sp)] = v
        species_n[escape_newick_label(sp)] = len(vals)

    print(f"🧮 Aggregated metrics for {len(species_val)} species")
    return species_val, species_n

# -------- write GraPhlAn annotation file --------
def write_graphlan_annotations(out_path: str, species_val: dict, mode: str, scale_height: bool):
    if not species_val:
        raise SystemExit("No species values to annotate.")

    # Value range for color scaling
    vmin, vmax = 0.0, 1.0

    lines = []

    # Global aesthetics
    lines.append(f"title\t{TITLE}")
    lines.append("title_font_size\t16")
    lines.append("annotation_background_width\t0.08")
    lines.append("annotation_background_alpha\t0.15")
    lines.append("annotation_legend_font_size\t8")
    lines.append("class_legend_font_size\t8")
    lines.append("branch_thickness\t0.7")
    lines.append("branch_color\t#666666")
    lines.append("total_plotted_degrees\t340")  # leave some space for ring label

    # Ring meta
    ring_label = RING_LABEL_PERCENT if mode == "percent" else RING_LABEL_MEAN
    lines.append(f"ring_label\t{RING_LEVEL}\t{ring_label}")
    lines.append(f"ring_label_font_size\t{RING_LEVEL}\t11")
    lines.append("ring_internal_separator_thickness\t1.0")
    lines.append("ring_external_separator_thickness\t1.0")
    lines.append("ring_separator_color\t#DDDDDD")

    # Base ring width/height
    # width: how thick the ring band is for each leaf (fraction of ring lane)
    lines.append(f"ring_width\t{RING_LEVEL}\t1.0")
    # default uniform height; override per-leaf below if SCALE_RING_HEIGHT
    lines.append(f"ring_height\t{RING_LEVEL}\t0.12")

    # Annotate each leaf: ring_color (gradient) + optional height + node marker color
    for leaf, v in species_val.items():
        color = value_to_hex(v, vmin, vmax)
        # ring color per leaf
        lines.append(f"{leaf}\tring_color\t{RING_LEVEL}\t{color}")
        if scale_height:
            # heights typically ~0.02–0.2; scale softly so 0 isn't invisible
            h = 0.06 + 0.14 * max(0.0, min(1.0, v))
            lines.append(f"{leaf}\tring_height\t{RING_LEVEL}\t{h}")
        # also color the leaf marker
        lines.append(f"{leaf}\tclade_marker_color\t{color}")
        lines.append(f"{leaf}\tclade_marker_size\t24")

    Path(out_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"📝 Wrote GraPhlAn annotations → {out_path}")

def main():
    ap = argparse.ArgumentParser(description="Build GraPhlAn annotations from prediction_files RNA scores.")
    ap.add_argument("--pred-dir", default=PREDICTION_DIR)
    ap.add_argument("--tree", default=TREE_NEWICK)
    ap.add_argument("--out-annot", default=ANNOT_OUT)
    ap.add_argument("--mode", choices=["percent","mean"], default=MODE)
    ap.add_argument("--threshold", type=float, default=THRESH)
    ap.add_argument("--scale-height", action="store_true", default=SCALE_RING_HEIGHT)
    args = ap.parse_args()

    # Quick sanity check that tree exists (so labels/escaping match your pipeline)
    if not Path(args.tree).exists():
        print(f"⚠️ Tree not found at {args.tree}. I’ll still write annotations, but make sure labels match your .nwk.")
    species_val, species_n = collect_species_values(args.pred_dir, args.mode, args.threshold)
    write_graphlan_annotations(args.out_annot, species_val, args.mode, args.scale_height)

    print("\nNext steps to render with GraPhlAn:")
    print(f"  graphlan_annotate.py {args.tree} tree_annot.xml --annot {args.out_annot}")
    print(f"  graphlan.py tree_annot.xml tree.png --dpi 300 --size 10")
    print("\nTips:")
    print("  • If some leaves don’t color, their labels in the .nwk likely differ (e.g., spacing).")
    print("    Open the .nwk and check a few leaf names; they should be your Organism with spaces → underscores.")
    print("  • Switch between modes with --mode percent|mean (and adjust --threshold for percent mode).")
    print("  • Add --scale-height if you want bar height to encode value as well.")
    print("  • You can tweak colors by editing C0/C1 in the script.")

if __name__ == "__main__":
    main()
