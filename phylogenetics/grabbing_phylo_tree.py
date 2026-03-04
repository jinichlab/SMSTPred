import re
from io import StringIO
from pathlib import Path
from collections import Counter
from typing import Dict, List, Optional, Tuple
import tempfile
import subprocess
import pandas as pd
from Bio import Phylo
import csv
import shlex

# === CONFIG ===
DROP_SUPERROOTS = {"cellular organisms"}
USE_SPECIES_AS_LEAVES = True
BRANCH_LENGTH = 1

CANONICAL_RANKS = [
    "domain", "kingdom", "phylum", "class",
    "order", "family", "genus", "species"
]

SAFE_CHARS_RE = re.compile(r"[^A-Za-z0-9_.|\-]+")

# === HELPERS ===
def split_lineage(s: str) -> List[str]:
    """Split a lineage string by comma or semicolon."""
    if not isinstance(s, str) or not s.strip():
        return []
    sep = "," if s.count(",") >= s.count(";") else ";"
    return [p.strip() for p in s.split(sep) if p.strip()]

def clean_taxon_name(tag: str) -> str:
    """Remove rank parentheses and clean spaces/quotes."""
    if not isinstance(tag, str):
        return ""
    m = re.match(r"^(.*?)(?:\s*\([^)]*\))?$", tag.strip())
    name = m.group(1) if m else tag
    return re.sub(r"\s+", " ", name.replace('"', "").replace("'", "")).strip()

def lineage_to_nodes(lineage_str: str) -> List[str]:
    """Return canonical taxon names from a lineage string."""
    nodes = []
    for raw in split_lineage(lineage_str):
        m = re.match(r"^(.*?)(?:\s*\(([^)]+)\))?$", raw)
        if not m:
            continue
        name, rank = clean_taxon_name(m.group(1)), (m.group(2) or "").lower()
        if name and rank in CANONICAL_RANKS and rank != "species":
            nodes.append(name)
    return [n for n in nodes if n not in DROP_SUPERROOTS]

def lineage_score(lineage_str: str) -> int:
    """Count canonical ranks in lineage."""
    return len(lineage_to_nodes(lineage_str))

def normalize_label(name: str) -> str:
    """Make label safe for Newick/GraPhlAn."""
    if not name:
        return "unnamed"
    s = str(name).strip().strip("'\"").replace("\xa0", " ")
    s = SAFE_CHARS_RE.sub("_", s)
    return re.sub(r"_+", "_", s).strip("_") or "unnamed"

# === TREE CLASSES ===
class TreeNode:
    __slots__ = ("name", "children", "is_leaf_tag")
    def __init__(self, name=None):
        self.name = name
        self.children: Dict[str, "TreeNode"] = {}
        self.is_leaf_tag = False

    def add_path(self, path, leaf_label=None):
        node = self
        for n in path:
            node = node.children.setdefault(n, TreeNode(n))
        if leaf_label:
            leaf = node.children.setdefault(leaf_label, TreeNode(leaf_label))
            leaf.is_leaf_tag = True

def to_newick(node: TreeNode, is_root=False) -> str:
    """Render Newick with constant branch length."""
    if not node.children:
        return f"{normalize_label(node.name)}:{BRANCH_LENGTH}"
    subs = ",".join(to_newick(c) for c in node.children.values())
    label = normalize_label(node.name) if node.name else ""
    return f"({subs}){label}" if is_root else f"({subs}){label}:{BRANCH_LENGTH}"

def build_tree_from_pairs(pairs: List[Tuple[str, str]]) -> TreeNode:
    """Build TreeNode tree from (species, lineage) pairs."""
    root, dup = TreeNode(None), Counter()
    for sp, lin in pairs:
        path = lineage_to_nodes(lin)
        if not path:
            continue
        if USE_SPECIES_AS_LEAVES:
            base = normalize_label(sp)
            dup[base] += 1
            leaf = base if dup[base] == 1 else f"{base}__{dup[base]}"
            root.add_path(path, leaf_label=leaf)
        else:
            root.add_path(path)
    return root

# === MAIN PIPELINE ===
def build_phylo_from_predictions(
    pred_dir: Path,
    score1_path: Path,
    score4_path: Path,
    out_phyloxml: Path = Path("taxonomy_tree.xml"),
    out_newick: Optional[Path] = Path("taxonomy_tree.nwk"),
    out_used: Path = Path("used_lineages.csv"),
    out_no_score14: Path = Path("species_without_lineage.csv"),
    verbose: bool = False,
):
    """Build phylogenetic tree from predictions + score tables."""
    # Load score-1
    score1 = pd.read_csv(score1_path, sep="\t", compression="gzip", dtype=str)
    entry2lin1 = dict(zip(score1["Entry"], score1["Taxonomic lineage"]))

    # Load score-4
    score4 = pd.read_csv(score4_path, sep="\t", compression="gzip", dtype=str)
    entry2lin4 = dict(zip(score4["Entry"], score4["Taxonomic lineage"]))

    # Gather species → entries
    species_to_entries: Dict[str, set] = {}
    for p in pred_dir.glob("*.csv"):
        df = pd.read_csv(p, dtype=str, encoding_errors="ignore", on_bad_lines="skip")
        if not {"Organism", "Entry"}.issubset(df.columns):
            continue
        for org, group in df.groupby("Organism"):
            species_to_entries.setdefault(org, set()).update(group["Entry"])

    # Pick best lineage per species
    used, no_lineage = [], []
    for sp, entries in species_to_entries.items():
        candidates = []
        for e in entries:
            if e in entry2lin1:
                candidates.append(("score1", entry2lin1[e]))
            elif e in entry2lin4:
                candidates.append(("score4", entry2lin4[e]))
        if not candidates:
            no_lineage.append(sp)
            continue
        candidates.sort(key=lambda x: (lineage_score(x[1]), x[0] == "score1"), reverse=True)
        used.append((sp, candidates[0][1], candidates[0][0]))

    # Save audit
    pd.DataFrame(used, columns=["species", "lineage", "source"]).to_csv(out_used, index=False)
    pd.DataFrame({"species": no_lineage}).to_csv(out_no_score14, index=False)

    # Build tree
    tree_root = build_tree_from_pairs([(sp, lin) for sp, lin, _ in used])
    newick = to_newick(tree_root, is_root=True) + ";"

    if out_newick:
        out_newick.write_text(newick)

    tree = Phylo.read(StringIO(newick), "newick")
    Phylo.write(tree, out_phyloxml, "phyloxml")

    if verbose:
        depths = [len(tree.get_path(leaf)) for leaf in tree.get_terminals()]
        print(f"Leaves: {len(depths)}  |  Max depth: {max(depths)}")

    return {"used": used, "newick": newick, "phyloxml_path": out_phyloxml}

import re, subprocess, shlex, tempfile, csv
from pathlib import Path
from typing import List, Tuple, Dict
from Bio import Phylo

# === CONSTANTS ===
HEX6 = re.compile(r"^#?[0-9A-Fa-f]{6}$")

GLOBAL_KEYS = {
    "title","title_font_size","branch_thickness","branch_color",
    "total_plotted_degrees","start_rotation",
}
RING_GLOBAL_KEYS = {
    "ring_label","ring_label_font_size","ring_width","ring_height",
    "ring_internal_separator_thickness","ring_external_separator_thickness","ring_separator_color",
}

# === UTILS ===
def _load_xml_leafset(xml_path: Path) -> set:
    """Return set of leaf names from a PhyloXML tree."""
    tree = Phylo.read(str(xml_path), "phyloxml")
    return {cl.name for cl in tree.get_terminals() if cl.name}

def strip_to_species_label(label: str) -> str:
    """Reduce label to Genus + species, drop strain/isolate info."""
    if not isinstance(label, str):
        return ""
    s = label.replace("_", " ").strip()
    s = re.split(r"\b(strain|isolate|DSM|ATCC|NCIMB|JCM)\b", s, flags=re.IGNORECASE)[0]
    parts = s.split()
    return "_".join(parts[:2]) if len(parts) >= 2 else (parts[0] if parts else "unnamed")

def map_to_existing_leaf(name: str, xml_leaves: set) -> Optional[str]:
    """Map candidate name to an existing XML leaf (underscore/space tolerant)."""
    for candidate in (name, name.replace(" ", "_"), name.replace("_", " ")):
        if candidate in xml_leaves:
            return candidate
    return None

def _parse_annot_file(annot_path: Path, xml_leaves: set):
    """Parse annotation file, normalize/match leaf labels."""
    globs, ring_globs, marker_lines, ring_lines = [], [], [], []
    for i, raw in enumerate(annot_path.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split("\t")

        # Normalize leaf names where applicable
        if parts[0] not in {"title","ring_label","ring_width","ring_height"}:
            candidate = strip_to_species_label(parts[0])
            mapped = map_to_existing_leaf(candidate, xml_leaves)
            if not mapped:
                continue
            parts[0] = mapped

        line = "\t".join(parts)
        if parts[0] in GLOBAL_KEYS:
            globs.append(line)
        elif parts[0] in RING_GLOBAL_KEYS:
            ring_globs.append(line)
        elif len(parts) >= 2 and parts[1] == "ring_color":
            ring_lines.append((i, line))
        else:
            marker_lines.append((i, line))

    return globs, ring_globs, marker_lines, ring_lines

def _precheck_line(ln_no: int, ln: str, xml_leaves: set) -> Tuple[bool, str]:
    """Lightweight syntax check before applying annotation line."""
    parts = ln.split("\t")
    if not (3 <= len(parts) <= 4):
        return False, f"bad_field_count({len(parts)})"
    leaf, opt = parts[0], parts[1]
    if leaf not in xml_leaves:
        return False, "leaf_not_in_xml"
    if len(parts) == 3:
        if opt in ("clade_marker_color","ring_color") and not HEX6.match(parts[2]):
            return False, "not_hex_color"
    elif opt != "ring_color" or not parts[2].isdigit() or not HEX6.match(parts[3]):
        return False, "invalid_ring_line"
    return True, ""

def _run_annot_with_lines(in_xml: Path, lines: List[str], out_xml: Path) -> bool:
    """Try applying annotation lines with graphlan_annotate."""
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as tmp:
        tmp.write("\n".join(lines) + "\n")
        tmp_path = Path(tmp.name)
    cmd = f"graphlan_annotate.py --annot {tmp_path} {in_xml} {out_xml}"
    subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    return out_xml.exists() and out_xml.stat().st_size > 0

def _bisect_apply(xml_start: Path, items: List[Tuple[int,str]], out_xml: Path,
                  tag: str, skipped: List[Dict], xml_leaves: set, min_chunk: int=1) -> Path:
    """Bisect apply annotation lines to isolate failures."""
    if items and _run_annot_with_lines(xml_start, [ln for _,ln in items], out_xml):
        return out_xml
    if not items:
        xml_start.replace(out_xml)
        return out_xml
    if len(items) <= min_chunk:
        for ln_no, ln in items:
            ok, why = _precheck_line(ln_no, ln, xml_leaves)
            skipped.append({"kind": tag, "line_no": ln_no, "line": ln,
                            "reason": ("annotate_failed" if ok else why)})
        xml_start.replace(out_xml)
        return out_xml
    mid = len(items)//2
    left, right = items[:mid], items[mid:]
    tmpL, tmpR = out_xml.parent / f"{out_xml.stem}.L.xml", out_xml.parent / f"{out_xml.stem}.R.xml"
    xml_after_left = _bisect_apply(xml_start, left, tmpL, tag, skipped, xml_leaves, min_chunk)
    xml_after_right = _bisect_apply(xml_after_left, right, tmpR, tag, skipped, xml_leaves, min_chunk)
    xml_after_right.replace(out_xml)
    for tmp in (tmpL, tmpR):
        tmp.unlink(missing_ok=True)
    return out_xml

# === PUBLIC API ===
def apply_annotations_with_diagnostics(xml_in: Path, annot_path: Path,
                                       xml_out: Path, batch_size=500,
                                       report_path: Optional[Path] = None):
    """Apply annotations incrementally, logging skipped lines."""
    if report_path is None:
        report_path = xml_out.with_suffix(".skipped.tsv")

    xml_leaves = _load_xml_leafset(xml_in)
    globs, ring_globs, marker_items, ring_items = _parse_annot_file(annot_path, xml_leaves)

    # globals
    cur = xml_in
    for tag, lines in (("global", globs), ("ring_globals", ring_globs)):
        nxt = xml_out.parent / f"{xml_out.stem}.{tag}.xml"
        if not _run_annot_with_lines(cur, lines or ["#"], nxt):
            raise RuntimeError(f"Failed applying {tag}")
        cur = nxt

    # markers
    skipped: List[Dict] = []
    for i in range(0, len(marker_items), batch_size):
        chunk = marker_items[i:i+batch_size]
        cur = _bisect_apply(cur, chunk, xml_out.parent / f"{xml_out.stem}.m{i}.xml",
                            "marker", skipped, xml_leaves)

    # rings
    for i in range(0, len(ring_items), batch_size):
        chunk = ring_items[i:i+batch_size]
        cur = _bisect_apply(cur, chunk, xml_out.parent / f"{xml_out.stem}.r{i}.xml",
                            "ring", skipped, xml_leaves)

    cur.replace(xml_out)

    # skipped log
    with report_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["kind","line_no","leaf","option","ring_level","value","reason","raw_line"])
        for rec in skipped:
            parts = rec["line"].split("\t")
            leaf, opt, lvl, val = (parts+[None]*4)[:4]
            w.writerow([rec["kind"], rec["line_no"], leaf, opt, lvl, val, rec["reason"], rec["line"]])

    print(f"✅ Annotation complete → {xml_out}")
    print(f"🧾 Skipped lines report → {report_path} ({len(skipped)} rows)")

    # --- Cleanup intermediate files ---
    for tmp in xml_out.parent.glob(f"{xml_out.stem}.m*.xml"):
        if tmp.name != xml_out.name:
            tmp.unlink(missing_ok=True)
    for tmp in xml_out.parent.glob(f"{xml_out.stem}.r*.xml"):
        if tmp.name != xml_out.name:
            tmp.unlink(missing_ok=True)

def render_graphlan(xml_in: Path, fig_out: Path, dpi=300, size=7.5, pad=0.24):
    """Render annotated tree with GraPhlAn."""
    cmd = f"graphlan.py {xml_in} {fig_out} --dpi {dpi} --size {size} --pad {pad} --avoid_reordering"
    print(f"🎨 Rendering: {cmd}")
    subprocess.run(shlex.split(cmd), check=True)
    print(f"✅ Figure written: {fig_out}")

def check_lineage_depth(xml_path: Path, max_allowed=8):
    """Check depth of leaves in tree."""
    tree = Phylo.read(str(xml_path), "phyloxml")
    depths = {t.name: tree.distance(t) for t in tree.get_terminals() if t.name}
    offenders = [name for name, d in depths.items() if d > max_allowed]
    print(f"🌿 Total leaves: {len(depths)}")
    print(f"📏 Max lineage depth: {max(depths.values()) if depths else 0}")
    print(f"⚠️ {len(offenders)} leaves exceed {max_allowed} levels")
    return offenders

from pathlib import Path

result = build_phylo_from_predictions(
    pred_dir=Path("../prediction_files"),   # folder with your *.csv prediction files
    score1_path=Path("../Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03.tsv.gz"),
    score4_path=Path("../Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03 (1).tsv.gz"),
    out_phyloxml=Path("taxonomy_tree.xml"),
    out_newick=Path("taxonomy_tree.nwk"),
    out_used=Path("used_lineages.csv"),
    out_no_score14=Path("species_without_lineage.csv"),
    verbose=True
)

print("Newick tree written:", result["newick"][:200], "...")
print("PhyloXML written to:", result["phyloxml_path"])

import random, re

def generate_graphlan_annotation(used_lineages_csv: Path, out_path: Path):
    """Generate a genus-colored GraPhlAn annotation file."""
    df = pd.read_csv(used_lineages_csv)

    # extract genus names
    genera = []
    for lin in df["lineage"].dropna():
        parts = [p.strip() for p in lin.split(";") if p.strip()]
        for p in reversed(parts):
            if "(genus)" in p.lower():
                genus = re.sub(r"\s*\([^)]*\)", "", p).strip()
                genera.append(genus)
                break
    genera = sorted(set(genera))

    random.seed(42)
    colors = ["#%06x" % random.randint(0, 0xFFFFFF) for _ in genera]

    with open(out_path, "w") as f:
        f.write("# Auto-generated GraPhlAn annotation file\n")
        f.write("title\tTaxonomy Tree\n")
        f.write("title_font_size\t20\n")
        f.write("branch_thickness\t0.4\n")
        f.write("total_plotted_degrees\t360\n")
        f.write("start_rotation\t90\n\n")
        for genus, color in zip(genera, colors):
            f.write(f"{genus}\tclade_marker_color\t{color}\n")
            f.write(f"{genus}\tclade_marker_size\t10\n")
    print(f"✅ Auto-generated annotation file with {len(genera)} genera: {out_path}")

import os
import re
import pandas as pd
from pathlib import Path

def escape_newick_label(label: str) -> str:
    if label is None:
        return ""
    label = label.replace(" ", "_")
    label = re.sub(r"[,:;()]", "_", label)
    return label

def lerp(a, b, t): return int(round(a + (b - a) * t))
def hex_to_rgb(hexcolor: str): return tuple(int(hexcolor.lstrip("#")[i:i+2], 16) for i in (0, 2, 4))
def rgb_to_hex(rgb): return "#{:02X}{:02X}{:02X}".format(*rgb)
C0, C1 = "#e9eef3", "#084594"

def value_to_hex(v: float, vmin=0.0, vmax=1.0) -> str:
    if v is None:
        return "#BDBDBD"  # neutral gray for missing
    v = max(vmin, min(vmax, v))
    t = 0.0 if vmax == vmin else (v - vmin) / (vmax - vmin)
    r0, g0, b0 = hex_to_rgb(C0)
    r1, g1, b1 = hex_to_rgb(C1)
    rgb = (lerp(r0, r1, t), lerp(g0, g1, t), lerp(b0, b1, t))  # ✅ include blue
    return rgb_to_hex(rgb)

def collect_species_values(pred_dir: str, mode="percent", thresh=0.5):
    """Aggregate per-species RNA-binding values."""
    species_records = {}
    for fname in os.listdir(pred_dir):
        if not fname.endswith(".csv"):
            continue
        fpath = os.path.join(pred_dir, fname)
        try:
            df = pd.read_csv(fpath, dtype=str, encoding_errors="ignore", on_bad_lines="skip")
        except Exception:
            continue
        org_col = next((c for c in ["Organism","organism","Species","species"] if c in df.columns), None)
        rna_col = next((c for c in ["RNA","rna","Rna"] if c in df.columns), None)
        if org_col is None or rna_col is None:
            continue
        df = df[[org_col, rna_col]].dropna()
        df[rna_col] = pd.to_numeric(df[rna_col], errors="coerce")
        for org, group in df.groupby(org_col):
            vals = group[rna_col].dropna().tolist()
            if vals:
                species_records.setdefault(org, []).extend(vals)

    species_val = {}
    for sp, vals in species_records.items():
        if not vals: continue
        v = float(pd.Series(vals).mean()) if mode == "mean" else float((pd.Series(vals) >= thresh).mean())
        species_val[escape_newick_label(sp)] = v
    return species_val

def write_graphlan_annotations(species_val: dict, out_path: str, mode="percent", ring_level=1):
    """Write a GraPhlAn annotation file with per-leaf ring colors."""
    vmin, vmax = 0.0, 1.0
    lines = [
        "title\tRNA-binding across species",
        "branch_thickness\t0.7",
        "branch_color\t#666666",
        "total_plotted_degrees\t340",
        f"ring_label\t{ring_level}\tRNA-binding %",
        f"ring_label_font_size\t{ring_level}\t11",
        "ring_internal_separator_thickness\t1.0",
        "ring_external_separator_thickness\t1.0",
        "ring_separator_color\t#DDDDDD",
        f"ring_width\t{ring_level}\t1.0",
        f"ring_height\t{ring_level}\t0.12"
    ]
    for leaf, v in species_val.items():
        color = value_to_hex(v, vmin, vmax)
        lines.append(f"{leaf}\tring_color\t{ring_level}\t{color}")
        lines.append(f"{leaf}\tclade_marker_color\t{color}")
        lines.append(f"{leaf}\tclade_marker_size\t24")
    Path(out_path).write_text("\n".join(lines) + "\n")
    print(f"✅ Wrote RNA-based annotation file → {out_path}")

XML_IN = Path("taxonomy_tree.xml")
ANNOT_TXT = Path("graphlan_annotations.txt")   #your annotation file
XML_ANNOTATED = Path("taxonomy_tree_annotated.xml")

species_val = collect_species_values("../prediction_files", mode="percent", thresh=0.5)
write_graphlan_annotations(species_val, ANNOT_TXT)

apply_annotations_with_diagnostics(
    xml_in=XML_IN,
    annot_path=ANNOT_TXT,
    xml_out=XML_ANNOTATED,
    batch_size=500,   # safe chunk size
    report_path=Path("taxonomy_tree.skipped.tsv")  #logs problematic lines
)

FIG_OUT = Path("tree_circular.png")

render_graphlan(
    xml_in=XML_ANNOTATED,
    fig_out=FIG_OUT,
    dpi=300,
    size=7.5,
    pad=0.24
)

def build_filtered_tree(genus_term: str,
                        used_csv: Path,
                        pred_dir: Path,
                        out_prefix: str):
    """Rebuild a genus-specific tree with proper coloring."""
    df = pd.read_csv(used_csv)
    df_filt = df[df["lineage"].str.contains(genus_term, case=False, na=False)]
    if df_filt.empty:
        print(f"⚠️ No entries found for '{genus_term}'")
        return

    df_filt = df_filt.copy()  # ✅ avoid SettingWithCopyWarning
    print(f"✅ Found {len(df_filt)} species for genus '{genus_term}'")

    # ✅ Trim the lineage to start at the genus level
    trimmed = []
    for lin in df_filt["lineage"]:
        parts = [p.strip() for p in lin.split(";") if p.strip()]
        for i, p in enumerate(parts):
            if genus_term.lower() in p.lower():
                parts = parts[i:]  # drop everything above the genus
                break
        trimmed.append("; ".join(parts))
    df_filt.loc[:, "lineage"] = trimmed

    # ✅ Build the subtree directly
    pairs = [(sp, lin) for sp, lin in zip(df_filt["species"], df_filt["lineage"])]
    from grabbing_phylo_tree import build_tree_from_pairs, to_newick

    # 🔧 remove 'Bacteria' filter in build_tree_from_pairs or replace with:
    # if not path: continue
    tree_root = build_tree_from_pairs(pairs)
    newick = to_newick(tree_root, is_root=True) + ";"

    out_xml = Path(f"{out_prefix}_tree.xml")
    tree = Phylo.read(StringIO(newick), "newick")
    Phylo.write(tree, out_xml, "phyloxml")

    # ✅ Write new annotation file matching these species
    from grabbing_phylo_tree import collect_species_values, write_graphlan_annotations
    species_val = collect_species_values(pred_dir, mode="percent", thresh=0.5)
    annot_txt = Path(f"{out_prefix}_annotations.txt")
    write_graphlan_annotations(species_val, annot_txt)

    # ✅ Annotate (with error handling)
    annot_xml = Path(f"{out_prefix}_annotated.xml")
    fig_png = Path(f"{out_prefix}_tree.png")

    print(f"🪶 Annotating: graphlan_annotate.py --annot {annot_txt} {out_xml} {annot_xml}")
    proc = subprocess.run(
        shlex.split(f"graphlan_annotate.py --annot {annot_txt} {out_xml} {annot_xml}"),
        capture_output=True, text=True
    )
    if not annot_xml.exists():
        print("⚠️ Annotation failed, falling back to unannotated tree")
        print(proc.stderr or proc.stdout)
        annot_xml = out_xml

    print(f"🎨 Rendering: graphlan.py {annot_xml} {fig_png}")
    subprocess.run(
        shlex.split(f"graphlan.py {annot_xml} {fig_png} --dpi 300 --size 7 --pad 0.24 --avoid_reordering"),
        check=False
    )

    print(f"✅ Done! Figure → {fig_png}")
    return fig_png

# Example usage:
build_filtered_tree("Mycobacterium",
                    used_csv=Path("used_lineages.csv"),
                    pred_dir=Path("../prediction_files"),
                    out_prefix="Mycobacterium")