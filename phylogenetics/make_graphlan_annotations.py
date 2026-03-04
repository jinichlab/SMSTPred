# ======= CONFIG =======
from pathlib import Path
import shlex, subprocess, re
from Bio import Phylo
import pandas as pd

PREDICTION_DIR = Path("prediction_files")
XML_IN          = Path("taxonomy_tree.xml")
ANNOT_TXT       = Path("graphlan_annotations.txt")
XML_ANNOTATED   = Path("taxonomy_tree_annotated.xml")
FIG_OUT         = Path("tree_circular.png")   # or .svg

TITLE = "RNA-binding across species (SAM MTases)"
MODE   = "percent"   # or "mean"
THRESH = 0.5

RING_LEVEL  = 1
RING_HEIGHT = 0.12
RING_WIDTH  = 1.0
LABEL_FONT  = 11
TITLE_FONT  = 16
TOTAL_DEGREES = 340
START_ROT     = 90

# ======= HELPERS =======
SAFE_CHARS_RE = re.compile(r"[^A-Za-z0-9_.|\-]+")

def normalize_leaf_label(name: str) -> str:
    if name is None: return ""
    s = str(name).strip()
    if (s.startswith("'") and s.endswith("'")) or (s.startswith('"') and s.endswith('"')):
        s = s[1:-1]
    s = s.replace("\xa0", " ")
    s = SAFE_CHARS_RE.sub("_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "unnamed"

def load_xml_leaves(xml_path: Path):
    tree = Phylo.read(str(xml_path), "phyloxml")
    leaves = [t.name for t in tree.get_terminals() if t.name]
    return set(leaves), leaves

def expand_to_existing_leaves(species_val, xml_leaves_set):
    out, misses = {}, []
    for raw_name, v in species_val.items():
        base = normalize_leaf_label(raw_name)
        if base in xml_leaves_set:
            out[base] = v
            continue
        pat = re.compile(rf"^{re.escape(base)}__(\d+)$")
        hit = False
        for leaf in xml_leaves_set:
            if leaf == base or pat.match(leaf):
                out[leaf] = v
                hit = True
        if not hit:
            misses.append((raw_name, base))
    return out, misses

def lerp(a,b,t): return int(round(a + (b-a)*t))
def hex_to_rgb(h): 
    h=h.lstrip("#"); 
    return tuple(int(h[i:i+2],16) for i in (0,2,4))
def rgb_to_hex(rgb): return "#{:02X}{:02X}{:02X}".format(*rgb)

C0, C1 = "#e9eef3", "#084594"
def value_to_hex(v, vmin=0.0, vmax=1.0):
    if v is None: return "#BDBDBD"
    v = max(vmin, min(vmax, v))
    t = 0.0 if vmax==vmin else (v-vmin)/(vmax-vmin)
    r0,g0,b0 = hex_to_rgb(C0); r1,g1,b1 = hex_to_rgb(C1)
    return rgb_to_hex((lerp(r0,r1,t), lerp(g0,g1,t), lerp(b0,b1,t)))

# ======= METRIC AGGREGATION =======
from typing import Dict, Tuple

def collect_species_values(pred_dir: Path, mode: str, thresh: float) -> Tuple[Dict[str,float], Dict[str,int]]:
    species_records = {}
    files = 0
    for p in sorted(pred_dir.glob("*.csv")):
        try:
            df = pd.read_csv(p, dtype=str, encoding_errors="ignore", on_bad_lines="skip", low_memory=False)
        except Exception as e:
            print(f"⚠️ Could not read {p.name}: {e}")
            continue
        files += 1
        org_col = next((c for c in ["Organism","organism","Species","species"] if c in df.columns), None)
        rna_col = next((c for c in ["RNA","rna","Rna"] if c in df.columns), None)
        if org_col is None or rna_col is None:
            print(f"⚠️ Missing Organism/RNA in {p.name} (have: {list(df.columns)[:6]} …)")
            continue
        rna = pd.to_numeric(df[rna_col], errors="coerce")
        org = df[org_col].astype(str)
        sub = pd.DataFrame({"org": org, "rna": rna}).dropna()
        if sub.empty: continue
        for sp, grp in sub.groupby("org"):
            vals = grp["rna"].astype(float).tolist()
            if vals:
                species_records.setdefault(sp, []).extend(vals)
    print(f"📄 Scanned {files} CSV file(s) under {pred_dir}")

    species_val, species_n = {}, {}
    for sp, vals in species_records.items():
        if not vals: continue
        v = float(pd.Series(vals).mean()) if mode=="mean" else float((pd.Series(vals) >= thresh).mean())
        species_val[sp] = v  # keep RAW names; we normalize against XML leaves later
        species_n[sp]   = len(vals)
    print(f"🧮 Aggregated metrics for {len(species_val)} species")
    return species_val, species_n

# ======= WRITE ANNOTATIONS TARGETED TO XML LEAVES =======
def write_graphlan_annotations_for_xml(xml_path: Path, out_path: Path, species_val, mode, thresh,
                                       title=TITLE, ring_level=RING_LEVEL, ring_height=RING_HEIGHT,
                                       ring_width=RING_WIDTH, total_deg=TOTAL_DEGREES, start_rot=START_ROT,
                                       title_font=TITLE_FONT, label_font=LABEL_FONT):
    xml_leaves_set, _ = load_xml_leaves(xml_path)
    leaf_value_map, misses = expand_to_existing_leaves(species_val, xml_leaves_set)
    if misses:
        print(f"⚠️ {len(misses)} species names did not match any XML leaves (showing up to 10):")
        for r, b in misses[:10]:
            print("   •", r, "→", b)

    ring_label = f"RNA-binding % (≥{thresh:.2f})" if mode == "percent" else "Mean RNA score"

    lines = [
        f"title\t{title}",
        f"title_font_size\t{title_font}",
        "branch_thickness\t0.7",
        "branch_color\t#666666",
        f"total_plotted_degrees\t{total_deg}",
        f"start_rotation\t{start_rot}",
        f"ring_label\t{ring_level}\t{ring_label}",
        f"ring_label_font_size\t{ring_level}\t{label_font}",
        f"ring_width\t{ring_level}\t{ring_width}",
        f"ring_height\t{ring_level}\t{ring_height}",
    ]

    for leaf, v in leaf_value_map.items():
        color = value_to_hex(v, 0.0, 1.0)
        lines.append(f"{leaf}\tring_color\t{ring_level}\t{color}")
        lines.append(f"{leaf}\tclade_marker_color\t{color}")
        lines.append(f"{leaf}\tclade_marker_size\t24")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"📝 Wrote annotations → {out_path} (targets {len(leaf_value_map)} of {len(xml_leaves_set)} leaves)")

# ======= TOLERANT ANNOTATE + RENDER =======
def run_graphlan(xml_in: Path, annot_txt: Path, xml_out: Path, fig_out: Path,
                 dpi=300, size=7.5, pad=0.24, avoid_reordering=True):
    if xml_out.exists():
        xml_out.unlink()
    cmd_annot = f"graphlan_annotate.py --annot {shlex.quote(str(annot_txt))} {shlex.quote(str(xml_in))} {shlex.quote(str(xml_out))}"
    print("🔧 Annotating:", cmd_annot)
    res = subprocess.run(shlex.split(cmd_annot), capture_output=True, text=True)
    print("annotate return:", res.returncode)
    if res.stdout.strip(): print("annotate STDOUT:\n", res.stdout)
    if res.stderr.strip(): print("annotate STDERR:\n", res.stderr)

    if not xml_out.exists() or xml_out.stat().st_size == 0:
        raise RuntimeError("graphlan_annotate.py failed and no annotated XML was produced.")

    # optional: quick content sanity
    txt = xml_out.read_text(errors="ignore")
    print("Contains clade_marker_color?", "clade_marker_color" in txt,
          "Contains ring_color?", "ring_color" in txt or "<ring" in txt)

    extra = " --avoid_reordering" if avoid_reordering else ""
    cmd_render = f"graphlan.py {shlex.quote(str(xml_out))} {shlex.quote(str(fig_out))} --dpi {dpi} --size {size} --pad {pad}{extra}"
    print("🎨 Rendering:", cmd_render)
    subprocess.run(shlex.split(cmd_render), check=True)
    print(f"✅ Wrote figure → {fig_out}")

# ======= RUN PIPELINE =======
species_val, species_n = collect_species_values(PREDICTION_DIR, MODE, THRESH)

write_graphlan_annotations_for_xml(
    xml_path=XML_IN,
    out_path=ANNOT_TXT,
    species_val=species_val,
    mode=MODE,
    thresh=THRESH,
    title=TITLE,
    ring_level=RING_LEVEL,
    ring_height=RING_HEIGHT,
    ring_width=RING_WIDTH,
    total_deg=TOTAL_DEGREES,
    start_rot=START_ROT,
    title_font=TITLE_FONT,
    label_font=LABEL_FONT
)

# quick check counts
annot_lines = Path(ANNOT_TXT).read_text().splitlines()
n_ring = sum(1 for ln in annot_lines if ln.split("\t")[1:2] == ["ring_color"])
n_mark = sum(1 for ln in annot_lines if ln.split("\t")[1:2] == ["clade_marker_color"])
print("ring_color lines:", n_ring, "clade_marker_color lines:", n_mark)

run_graphlan(XML_IN, ANNOT_TXT, XML_ANNOTATED, FIG_OUT, dpi=300, size=7.5, pad=0.24, avoid_reordering=True)

# display
from IPython.display import Image, SVG, display
display(Image(filename=str(FIG_OUT)) if FIG_OUT.suffix.lower() != ".svg" else SVG(filename=str(FIG_OUT)))