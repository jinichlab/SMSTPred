import re
from io import StringIO
from pathlib import Path
from collections import Counter
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from Bio import Phylo

# === knobs you can tweak ===
DROP_SUPERROOTS = {"cellular organisms"}   # remove this super-root if present
USE_SPECIES_AS_LEAVES = True               # attach organism as leaf under last taxon node
BRANCH_LENGTH = 1                          # constant branch length for Newick

def split_lineage(s: str, sep_guess=None) -> List[str]:
    """Split a lineage string by comma or semicolon, robust to spaces."""
    if not isinstance(s, str):
        return []
    s = s.strip()
    if not s:
        return []
    sep = sep_guess
    if sep is None:
        c_comma = s.count(",")
        c_semi  = s.count(";")
        sep = "," if c_comma >= c_semi else ";"
    parts = [p.strip() for p in s.split(sep)]
    return [p for p in parts if p]

def strip_rank(tag: str) -> str:
    """Remove trailing ' (rank)' if present, keep the taxon name only."""
    if not isinstance(tag, str):
        return ""
    m = re.match(r"^(.*?)(?:\s*\([^)]*\))\s*$", tag)
    name = m.group(1) if m else tag
    name = re.sub(r"\s+", " ", name).strip()
    return name

def clean_taxon_name(name: str) -> str:
    """Basic cleaning of a taxon node label."""
    name = strip_rank(name)
    name = name.replace('"', "").replace("'", "").replace("\xa0", " ").strip()
    return name

def lineage_to_nodes(lineage_str: str) -> List[str]:
    """Turn a lineage string into a list of cleaned taxon names."""
    parts = split_lineage(lineage_str)
    parts = [clean_taxon_name(p) for p in parts]
    parts = [p for p in parts if p]
    # drop super roots (optional)
    parts = [p for p in parts if p not in DROP_SUPERROOTS]
    # de-duplicate accidental repeats like "...; Bacilli; Bacillota; Bacilli"
    dedup = []
    for p in parts:
        if not dedup or dedup[-1] != p:
            dedup.append(p)
    return dedup

class TreeNode:
    __slots__ = ("name", "children", "is_leaf_tag")
    def __init__(self, name=None):
        self.name = name
        self.children = {}  # name -> TreeNode
        self.is_leaf_tag = False

    def add_path(self, path, leaf_label=None):
        node = self
        for n in path:
            if n not in node.children:
                node.children[n] = TreeNode(n)
            node = node.children[n]
        if leaf_label is not None:
            if leaf_label not in node.children:
                node.children[leaf_label] = TreeNode(leaf_label)
            node.children[leaf_label].is_leaf_tag = True

SAFE_CHARS_RE = re.compile(r"[^A-Za-z0-9_.|\-]+")

def normalize_leaf_label(name: str) -> str:
    if name is None:
        return ""
    s = str(name).strip()
    # strip surrounding quotes
    if (s.startswith("'") and s.endswith("'")) or (s.startswith('"') and s.endswith('"')):
        s = s[1:-1]
    s = s.replace("\xa0", " ")
    s = SAFE_CHARS_RE.sub("_", s)      # turn weird chars (/, ', etc.) into _
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "unnamed"

def escape_newick_label(label: str) -> str:
    # Use the same normalization everywhere
    return normalize_leaf_label(label)

def to_newick(node: TreeNode, is_root=False) -> str:
    """Recursively render Newick with constant branch lengths; omit branch length on root."""
    if not node.children:
        label = escape_newick_label(node.name) if node.name is not None else ""
        return f"{label}:{BRANCH_LENGTH}"
    subs = ",".join(to_newick(child, is_root=False) for child in node.children.values())
    label = "" if node.name is None else escape_newick_label(node.name)
    if is_root:
        return f"({subs}){label}"
    else:
        return f"({subs}){label}:{BRANCH_LENGTH}"

def escape_newick_label_for_check(label: str) -> str:
    return normalize_leaf_label(label)

def build_phylo_from_predictions(
    pred_dir: Path,
    score1_path: Path,
    score4_path: Path,
    tree_check: Optional[Path] = None,
    out_used: Path = Path("used_lineages_from_score14.csv"),
    out_no_score14: Path = Path("species_without_score14_protein.csv"),
    out_not_in_tree: Path = Path("species_not_in_tree.csv"),
    out_phyloxml: Path = Path("taxonomy_tree.xml"),
    out_newick: Optional[Path] = None,  # set a Path if you also want a .nwk saved
    verbose: bool = False,
):
    """
    Jupyter-friendly version of your CLI:
      - scans prediction CSVs to map species -> Entry IDs
      - pulls lineage from score-1 (preferred) or score-4 (fallback)
      - builds a tree (constant branch length)
      - writes PhyloXML (and optionally Newick)
      - writes audit CSVs
    Returns a dict with summary + DataFrames.
    """
    pred_dir = Path(pred_dir)
    if not pred_dir.exists():
        raise FileNotFoundError(f"Prediction folder not found: {pred_dir}")

    # --- Load score-1 table ---
    print(f"📥 Loading score-1 table: {score1_path}")
    score1 = pd.read_csv(score1_path, sep="\t", compression="gzip", dtype=str, low_memory=False)
    needed1 = {"Entry", "Taxonomic lineage"}
    if not needed1.issubset(score1.columns):
        raise ValueError(f"{score1_path} must contain columns: {needed1} (found {list(score1.columns)[:10]} ...)")
    score1 = score1[["Entry", "Taxonomic lineage"]].dropna()
    score1["Entry"] = score1["Entry"].astype(str).str.strip()
    entry2lin1: Dict[str, str] = dict(zip(score1["Entry"], score1["Taxonomic lineage"]))
    print(f"✅ Loaded {len(entry2lin1):,} score-1 entries with lineage")

    # --- Load score-4 table (fallback) ---
    print(f"📥 Loading score-4 table: {score4_path}")
    score4 = pd.read_csv(score4_path, sep="\t", compression="gzip", dtype=str, low_memory=False)
    needed4 = {"Entry", "Taxonomic lineage"}
    if not needed4.issubset(score4.columns):
        raise ValueError(f"{score4_path} must contain columns: {needed4} (found {list(score4.columns)[:10]} ...)")
    score4 = score4[["Entry", "Taxonomic lineage"]].dropna()
    score4["Entry"] = score4["Entry"].astype(str).str.strip()
    entry2lin4: Dict[str, str] = dict(zip(score4["Entry"], score4["Taxonomic lineage"]))
    print(f"✅ Loaded {len(entry2lin4):,} score-4 entries with lineage")

    # --- Scan prediction files and collect per-species Entry IDs ---
    species_to_entries: Dict[str, set] = {}
    csv_files = sorted([p for p in pred_dir.iterdir() if p.suffix.lower() == ".csv"])
    if not csv_files:
        print(f"⚠️ No CSVs found under {pred_dir}")
    print(f"🔎 Scanning {len(csv_files)} prediction CSVs...")
    for p in csv_files:
        try:
            df = pd.read_csv(p, dtype=str, encoding_errors="ignore", on_bad_lines="skip", low_memory=False)
        except Exception as e:
            print(f"  ⚠️ Could not read {p.name}: {e}")
            continue
        if not {"Organism", "Entry"}.issubset(df.columns):
            if verbose:
                print(f"  ⚠️ Missing Organism/Entry in {p.name} (have: {list(df.columns)[:8]} …)")
            continue
        sub = df[["Organism", "Entry"]].dropna()
        sub["Organism"] = sub["Organism"].astype(str).str.strip()
        sub["Entry"] = sub["Entry"].astype(str).str.strip()
        for org, group in sub.groupby("Organism"):
            species_to_entries.setdefault(org, set()).update(group["Entry"].tolist())

    unique_species = len(species_to_entries)
    total_entries_seen = sum(len(v) for v in species_to_entries.values())
    print(f"🧾 Collected entries for {unique_species:,} species "
          f"({total_entries_seen:,} unique Entry IDs across all species)")

    # --- For each species: prefer score-1; fallback score-4; take first (sorted) lineage ---
    used_rows: List[Tuple[str, str, str, int, int]] = []  # (species, lineage, source, n1, n4)
    no_score14_species: List[str] = []
    matched_via_1 = 0
    matched_via_4 = 0

    for sp, entries in species_to_entries.items():
        matches1 = sorted([e for e in entries if e in entry2lin1])
        matches4 = sorted([e for e in entries if e in entry2lin4])
        if matches1:
            lineage = entry2lin1[matches1[0]]   # deterministic first
            used_rows.append((sp, lineage, "score1", len(matches1), len(matches4)))
            matched_via_1 += 1
        elif matches4:
            lineage = entry2lin4[matches4[0]]
            used_rows.append((sp, lineage, "score4", 0, len(matches4)))
            matched_via_4 += 1
        else:
            no_score14_species.append(sp)

    print(f"✅ Matched {matched_via_1 + matched_via_4:,} species total "
          f"({matched_via_1:,} via score-1, {matched_via_4:,} via score-4 fallback)")
    print(f"❌ Species with NO score-1/4 protein: {len(no_score14_species):,}")

    # --- (Optional) compare matched species to an existing tree (for reporting only) ---
    not_in_tree: List[str] = []
    if tree_check:
        tree_path = Path(tree_check)
        if not tree_path.exists():
            print(f"  ⚠️ Tree not found at {tree_path} — skipping tree comparison")
        else:
            nwk_text = tree_path.read_text(encoding="utf-8", errors="ignore")
            leaf_labels = set(re.findall(r"([^\s,():]+):", nwk_text))
            for sp, *_rest in used_rows:
                if escape_newick_label_for_check(sp) not in leaf_labels:
                    not_in_tree.append(sp)
            print(f"🌳 Matched species NOT present as leaves in provided tree: {len(not_in_tree):,}")

    # --- Write audits ---
    used_df = None
    if used_rows:
        used_df = pd.DataFrame(
            used_rows,
            columns=["species", "lineage", "source", "n_entries_matched_1", "n_entries_matched_4"]
        )
        used_df.to_csv(out_used, index=False)
        print(f"📝 Wrote used lineages → {out_used}")

    if no_score14_species:
        pd.DataFrame({"species": sorted(no_score14_species)}).to_csv(out_no_score14, index=False)
        print(f"🧭 Wrote species without score-1/4 protein → {out_no_score14}")

    if tree_check and not_in_tree:
        pd.DataFrame({"species": sorted(not_in_tree)}).to_csv(out_not_in_tree, index=False)
        print(f"🌱 Wrote matched-but-not-in-tree species → {out_not_in_tree}")

    # --- Build tree (your same logic) -> Newick string ---
    if not used_rows:
        print("❌ No matched organisms; tree not created.")
        return {
            "used_df": used_df,
            "no_score14_species": no_score14_species,
            "not_in_tree": not_in_tree,
            "newick": None,
            "phyloxml_path": None,
        }

    # Assemble pairs and build the internal tree
    used_pairs = [(sp, lineage) for (sp, lineage, *_rest) in used_rows]

    # Build TreeNode structure
    root = TreeNode(name=None)
    dup_counter = Counter()
    attached = 0
    for sp, lineage in used_pairs:
        path = lineage_to_nodes(lineage)
        if not path:
            continue
        if USE_SPECIES_AS_LEAVES:
            base = normalize_leaf_label(sp)
            dup_counter[base] += 1
            leaf = base if dup_counter[base] == 1 else f"{base}__{dup_counter[base]}"
            root.add_path(path, leaf_label=leaf)
        else:
            root.add_path(path)
        attached += 1
    if verbose:
        print(f"🌿 Attached {attached} leaves to tree")

    newick = to_newick(root, is_root=True)
    if not newick.endswith(";"):
        newick += ";"

    # --- Write optional .nwk ---
    if out_newick is not None:
        Path(out_newick).write_text(newick, encoding="utf-8")
        print(f"🌳 Wrote Newick → {out_newick}")

    # --- Convert Newick → PhyloXML and write ---
    tree = Phylo.read(StringIO(newick), "newick")
    Phylo.write(tree, out_phyloxml, "phyloxml")
    print(f"🧬 Wrote PhyloXML → {out_phyloxml}")

    return {
        "used_df": used_df,
        "no_score14_species": no_score14_species,
        "not_in_tree": not_in_tree,
        "newick": newick,
        "phyloxml_path": Path(out_phyloxml),
    }

def phyloxml_from_used_csv(
    used_csv: Path = Path("used_lineages_from_score14.csv"),
    out_phyloxml: Path = Path("taxonomy_tree.xml"),
    out_newick: Optional[Path] = None,
):
    """
    Jupyter refactor of your second script: reads the 'used' CSV
    and writes PhyloXML (and optionally Newick).
    """
    path = Path(used_csv)
    if not path.exists():
        raise FileNotFoundError(f"Used-lineages CSV not found: {used_csv}")

    df = pd.read_csv(path, dtype=str)
    cols = {c.lower(): c for c in df.columns}
    sp_col = cols.get("species") or cols.get("organism")
    lin_col = cols.get("lineage")
    if not sp_col or not lin_col:
        raise ValueError(f"{used_csv} must have columns including 'species' (or 'organism') and 'lineage'")

    pairs: List[Tuple[str, str]] = []
    for _, row in df.iterrows():
        sp = (row.get(sp_col) or "").strip()
        lin = (row.get(lin_col) or "").strip()
        if sp and lin:
            pairs.append((sp, lin))
    if not pairs:
        raise ValueError("No (species, lineage) pairs found in used CSV")

    # Build TreeNode structure
    root = TreeNode(name=None)
    dup = Counter()
    attached = 0
    for sp, lin in pairs:
        path = lineage_to_nodes(lin)
        if not path:
            continue
        if USE_SPECIES_AS_LEAVES:
            leaf = sp
            dup[sp] += 1
            if dup[sp] > 1:
                leaf = f"{sp}__{dup[sp]}"
            root.add_path(path, leaf_label=leaf)
        else:
            root.add_path(path)
        attached += 1

    newick = to_newick(root, is_root=True)
    if not newick.endswith(";"):
        newick += ";"

    if out_newick is not None:
        Path(out_newick).write_text(newick, encoding="utf-8")
        print(f"🌳 Wrote Newick → {out_newick} (attached {attached} leaves)")

    tree = Phylo.read(StringIO(newick), "newick")
    Phylo.write(tree, out_phyloxml, "phyloxml")
    print(f"🧬 Wrote PhyloXML → {out_phyloxml} (attached {attached} leaves)")

result = build_phylo_from_predictions(
    pred_dir=Path("prediction_files"),
    score1_path=Path("Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03.tsv.gz"),
    score4_path=Path("Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03 (1).tsv.gz"),
    tree_check=None,  # or Path("some_existing_tree.nwk")
    out_used=Path("used_lineages_from_score14.csv"),
    out_no_score14=Path("species_without_score14_protein.csv"),
    out_not_in_tree=Path("species_not_in_tree.csv"),
    out_phyloxml=Path("taxonomy_tree.xml"),
    out_newick=Path("taxonomy_tree.nwk"),  # optional
    verbose=True,
)

import subprocess, shlex, tempfile, os
from pathlib import Path
from IPython.display import Image, display, SVG

xml_path = "taxonomy_tree.xml"  # your phyloXML from earlier
out_png = "graphlan_tree.png"   # or .svg

# Basic render: GraPhlAn can render a plain phyloXML directly
cmd = f"graphlan.py {shlex.quote(xml_path)} {shlex.quote(out_png)} --dpi 300 --size 7 --pad 0.2"
print("Running:", cmd)
subprocess.run(shlex.split(cmd), check=True)

# Display in notebook
display(Image(filename=out_png))

# ======= IMPORTS & PATHS =======
from pathlib import Path
import shlex, subprocess, tempfile, re
from Bio import Phylo
import pandas as pd

PREDICTION_DIR = Path("prediction_files")
XML_IN         = Path("taxonomy_tree.xml")            # from your build step
ANNOT_TXT      = Path("graphlan_annotations.txt")     # written by your writer
XML_ANNOTATED  = Path("taxonomy_tree_annotated.xml")  # will be produced here
FIG_OUT        = Path("tree_circular.png")            # or .svg

# ======= INCREMENTAL ANNOTATION UTILS =======
# --- put this with your other helpers (after imports for the annotation step) ---
def _probe_graphlan_support() -> dict:
    """Return a dict of supported directives by trying a tiny fake tree."""
    import tempfile, shlex, subprocess
    from pathlib import Path

    # minimal PhyloXML with two leaves
    mini_xml = """<?xml version="1.0" encoding="UTF-8"?>
<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
          xmlns="http://www.phyloxml.org" xsi:schemaLocation="http://www.phyloxml.org
          http://www.phyloxml.org/1.10/phyloxml.xsd">
  <phylogeny rooted="false">
    <clade><clade><name>a</name></clade><clade><name>b</name></clade></clade>
  </phylogeny>
</phyloxml>"""
    tmp = Path(tempfile.mkdtemp())
    x = tmp / "t.xml"; x.write_text(mini_xml)
    results = {}

    def try_ann(lines):
        ann = tmp / "a.txt"; ann.write_text("\n".join(lines) + "\n")
        out = tmp / "o.xml"
        if out.exists(): out.unlink()
        cmd = f"graphlan_annotate.py --annot {shlex.quote(str(ann))} {shlex.quote(str(x))} {shlex.quote(str(out))}"
        r = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
        ok = out.exists() and out.stat().st_size > 0
        return ok, r.stderr.strip(), r.stdout.strip()

    # test each directive we might use
    tests = {
        "clade_marker_color": ["a\tclade_marker_color\t#FF0000"],
        "clade_marker_size":  ["a\tclade_marker_size\t24"],
        "ring_color":         ["a\tring_color\t1\t#FF0000",  # the problematic one on your machine
                               "ring_label\t1\tTest", "ring_width\t1\t1", "ring_height\t1\t0.1"],
    }

    for key, lines in tests.items():
        ok, err, out = try_ann(lines)
        results[key] = {"ok": ok, "stderr": err, "stdout": out}

    # clean up tempdir silently
    try:
        for p in tmp.iterdir(): p.unlink()
        tmp.rmdir()
    except Exception:
        pass

    return results

# ======= DIAGNOSTIC INCREMENTAL ANNOTATION =======
import re, tempfile, shlex, subprocess
from pathlib import Path
from typing import List, Tuple, Dict
from Bio import Phylo

HEX6 = re.compile(r"^#?[0-9A-Fa-f]{6}$")

GLOBAL_KEYS = {
    "title","title_font_size","branch_thickness","branch_color",
    "total_plotted_degrees","start_rotation",
}
RING_GLOBAL_KEYS = {
    "ring_label","ring_label_font_size","ring_width","ring_height",
    "ring_internal_separator_thickness","ring_external_separator_thickness","ring_separator_color",
}

def _load_xml_leafset(xml_path: Path) -> set:
    t = Phylo.read(str(xml_path), "phyloxml")
    return {cl.name for cl in t.get_terminals() if cl.name}

def _parse_annot_file(annot_path: Path) -> Tuple[List[str], List[str], List[Tuple[int,str]], List[Tuple[int,str]]]:
    """Return (global_lines, ring_global_lines, marker_lines[(ln_no,line)], ring_lines[(ln_no,line)])"""
    globs, ring_globs, marker_lines, ring_lines = [], [], [], []
    for i, raw in enumerate(annot_path.read_text(encoding="utf-8").splitlines(), start=1):
        ln = raw.strip("\n")
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if not parts:
            continue
        if parts[0] in GLOBAL_KEYS:
            globs.append(ln); continue
        if parts[0] in RING_GLOBAL_KEYS:
            ring_globs.append(ln); continue
        # per-leaf
        if len(parts) >= 2 and parts[1] == "ring_color":
            ring_lines.append((i, ln))
        else:
            marker_lines.append((i, ln))
    return globs, ring_globs, marker_lines, ring_lines

def _precheck_line(ln_no: int, ln: str, xml_leaves: set) -> Tuple[bool, str]:
    """
    Validate syntax & obvious issues before calling GraPhlAn.
    Returns (ok, reason_if_not_ok).
    """
    parts = ln.split("\t")
    n = len(parts)
    if n not in (3,4):
        return False, f"bad_field_count({n})"
    leaf = parts[0]
    if leaf not in xml_leaves:
        return False, "leaf_not_in_xml"
    opt = parts[1]
    if n == 3:
        # marker line: leaf  opt  value
        val = parts[2]
        if opt in ("clade_marker_color","ring_color"):  # ring_color shouldn't be here but just in case
            if not HEX6.match(val): return False, "not_hex_color"
    else:
        # ring line: leaf  ring_color  level  #RRGGBB
        level, color = parts[2], parts[3]
        if opt != "ring_color":
            return False, "unexpected_opt_for_ring"
        if not level.isdigit():
            return False, "ring_level_not_int"
        if not HEX6.match(color):
            return False, "not_hex_color"
    return True, ""

def _run_annot_with_lines(in_xml: Path, lines: List[str], out_xml: Path) -> bool:
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as tmp:
        tmp.write("\n".join(lines) + "\n")
        tmp_path = Path(tmp.name)
    cmd = f"graphlan_annotate.py --annot {shlex.quote(str(tmp_path))} {shlex.quote(str(in_xml))} {shlex.quote(str(out_xml))}"
    res = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    return out_xml.exists() and out_xml.stat().st_size > 0

def _bisect_apply(
    xml_start: Path,
    items: List[Tuple[int,str]],     # [(line_no, line), ...]
    out_xml: Path,
    tag: str,
    skipped: List[Dict],
    xml_leaves: set,
    min_chunk: int = 1,
) -> Path:
    """
    Try apply all items; if failing, bisect down to single lines,
    record those lines to `skipped` with a reason.
    Returns the path of an XML that includes all applied (non-toxic) changes.
    """
    # Try all at once
    if items:
        if _run_annot_with_lines(xml_start, [ln for _,ln in items], out_xml):
            return out_xml

    # Base cases
    if not items:
        # Nothing to apply; just pass through xml_start
        if out_xml.exists(): out_xml.unlink()
        xml_start.replace(out_xml)
        return out_xml

    if len(items) <= min_chunk:
        # If it's a single line (or small chunk), do precheck and then declare it toxic
        for ln_no, ln in items:
            ok, why = _precheck_line(ln_no, ln, xml_leaves)
            skipped.append({
                "kind": tag,
                "line_no": ln_no,
                "line": ln,
                "reason": ("annotate_failed" if ok else why),
            })
        # No change; copy input forward
        if out_xml.exists(): out_xml.unlink()
        xml_start.replace(out_xml)
        return out_xml

    # Split and recurse
    mid = len(items)//2
    left, right = items[:mid], items[mid:]
    tmpL = out_xml.parent / f"{out_xml.stem}.__{tag}_L__.xml"
    tmpR = out_xml.parent / f"{out_xml.stem}.__{tag}_R__.xml"

    xml_after_left = _bisect_apply(xml_start, left, tmpL, tag, skipped, xml_leaves, min_chunk)
    xml_after_right = _bisect_apply(xml_after_left, right, tmpR, tag, skipped, xml_leaves, min_chunk)

    # Move final into out_xml
    if out_xml.exists(): out_xml.unlink()
    xml_after_right.replace(out_xml)
    # cleanup leftovers best-effort
    for p in (tmpL, tmpR):
        try: p.unlink()
        except: pass
    return out_xml

def apply_annotations_with_diagnostics(xml_in: Path, annot_path: Path, xml_out: Path, batch_size=500, report_path: Path = None):
    """
    Applies annotations in a few stages and records every skipped (toxic) line to report_path.
    """
    if report_path is None:
        report_path = xml_out.with_suffix(".skipped.tsv")

    xml_leaves = _load_xml_leafset(xml_in)
    globs, ring_globs, marker_items, ring_items = _parse_annot_file(annot_path)

    # 1) globals -> xml_g
    xml_g = xml_out.parent / f"{xml_out.stem}.g.xml"
    if not _run_annot_with_lines(xml_in, globs or ["#"], xml_g):
        raise RuntimeError("Failed applying global options")

    # 2) ring globals -> xml_rg
    xml_rg = xml_out.parent / f"{xml_out.stem}.rg.xml"
    if not _run_annot_with_lines(xml_g, ring_globs or ["#"], xml_rg):
        raise RuntimeError("Failed applying ring global options")

    # 3) marker lines in batches with bisection
    skipped: List[Dict] = []
    cur = xml_rg
    print(f"Applying {len(marker_items)} marker lines in batches of {batch_size}…")
    for i in range(0, len(marker_items), batch_size):
        chunk = marker_items[i:i+batch_size]
        nxt = xml_out.parent / f"{xml_out.stem}.m{i//batch_size}.xml"
        cur = _bisect_apply(cur, chunk, nxt, tag="marker", skipped=skipped, xml_leaves=xml_leaves, min_chunk=1)

    # 4) ring_color lines in batches with bisection
    print(f"Applying {len(ring_items)} ring_color lines in batches of {batch_size}…")
    for i in range(0, len(ring_items), batch_size):
        chunk = ring_items[i:i+batch_size]
        nxt = xml_out.parent / f"{xml_out.stem}.r{i//batch_size}.xml"
        cur = _bisect_apply(cur, chunk, nxt, tag="ring", skipped=skipped, xml_leaves=xml_leaves, min_chunk=1)

    # finalize
    if xml_out.exists(): xml_out.unlink()
    cur.rename(xml_out)

    # write skipped report
    import csv
    with report_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["kind","line_no","leaf","option","ring_level","value","reason","raw_line"])
        for rec in skipped:
            parts = rec["line"].split("\t")
            leaf  = parts[0] if parts else ""
            opt   = parts[1] if len(parts) >= 2 else ""
            lvl   = parts[2] if len(parts) >= 3 else ""
            val   = parts[3] if len(parts) >= 4 else (parts[2] if len(parts) >= 3 else "")
            w.writerow([rec["kind"], rec["line_no"], leaf, opt, lvl, val, rec["reason"], rec["line"]])
    print(f"✅ Annotation complete → {xml_out}")
    print(f"🧾 Skipped lines report → {report_path} (rows: {sum(1 for _ in report_path.open())-1})")

apply_annotations_with_diagnostics(XML_IN, ANNOT_TXT, XML_ANNOTATED, batch_size=500)

cmd = f"graphlan.py {shlex.quote(str(XML_ANNOTATED))} {shlex.quote(str(FIG_OUT))} --dpi 300 --size 7.5 --pad 0.24 --avoid_reordering"
print("🎨 Rendering:", cmd)
subprocess.run(shlex.split(cmd), check=True)
print("✅ Figure written:", FIG_OUT)

# === SAFE XML PATCHER: promote failed ring_color to internal annotation rings ===
from pathlib import Path
import re, shlex, subprocess
import xml.etree.ElementTree as ET
import pandas as pd

XML_IN    = Path("taxonomy_tree_annotated.xml")          # produced by the diagnostic annotator
SKIPPED   = Path("taxonomy_tree_annotated.skipped.tsv")   # report it wrote
ANNOT_TXT = Path("graphlan_annotations.txt")              # original annotations
XML_OUT   = Path("taxonomy_tree_annotated.patched.xml")
PNG_OUT   = Path("tree_circular_patched.png")
BIG_PNG   = Path("tree_circular_patched_big.png")
SVG_OUT   = Path("tree_circular_patched.svg")

# --- load failed ring_color leaves ---
df = pd.read_csv(SKIPPED, sep="\t")
df = df[(df["kind"] == "ring")]
failed_leaves = df["leaf"].dropna().astype(str).tolist()
print(f"Found {len(failed_leaves)} ring failures")

if not failed_leaves:
    print("Nothing to patch; skipping.")
else:
    # --- build target color map from original annotation text ---
    target_colors = {}  # leaf -> color
    hex6 = re.compile(r"^#[0-9A-Fa-f]{6}$")
    for raw in ANNOT_TXT.read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#"): 
            continue
        parts = raw.split("\t")
        if len(parts) == 4 and parts[1] == "ring_color":
            leaf, _, level, color = parts
            if hex6.match(color):
                # keep first
                target_colors.setdefault(leaf, color)

    missing = [l for l in failed_leaves if l not in target_colors]
    if missing:
        print(f"⚠️ {len(missing)} failed leaves have no ring_color in annot file (will be skipped). e.g. {missing[:3]}")

    to_patch = [(l, target_colors[l]) for l in failed_leaves if l in target_colors]
    print(f"Will patch {len(to_patch)} leaves")

    # --- parse PhyloXML with namespaces ---
    tree = ET.parse(XML_IN)
    root = tree.getroot()

    # detect namespace (usual phyloxml ns)
    if root.tag.startswith("{"):
        ns_uri = root.tag.split("}", 1)[0][1:]
    else:
        ns_uri = "http://www.phyloxml.org"  # fallback
    NS = {"px": ns_uri}
    def q(tag):  # qualify tag with namespace
        return f"{{{ns_uri}}}{tag}"

    # a helper to find clade by its <name> text (exact match after your normalization)
    def find_clade_by_name(root, leaf_name):
        # search all <name> nodes and check their text
        for name_el in root.findall(".//px:name", NS):
            if (name_el.text or "") == leaf_name:
                # return the enclosing clade of this name node
                parent = name_el
                while parent is not None and parent.tag != q("clade"):
                    parent = parent.getparent() if hasattr(parent, "getparent") else None
                    # xml.etree doesn't keep parent pointers; walk up manually:
                    # fallback: do a small upward search by re-finding the clade that contains this name
                if parent is None:
                    # xml.etree lacks getparent; do a second pass to locate the clade
                    # by scanning clades and checking if their immediate name child is our node
                    for cl in root.findall(".//px:clade", NS):
                        nm = cl.find("px:name", NS)
                        if nm is not None and nm.text == leaf_name:
                            return cl
                return parent
        return None

    # Because xml.etree lacks .getparent(), we implement a direct clade lookup:
    def find_clade_by_name_fast(root, leaf_name):
        for cl in root.findall(".//px:clade", NS):
            nm = cl.find("px:name", NS)
            if nm is not None and (nm.text or "") == leaf_name:
                return cl
        return None

    # ensure an <annotation> child exists and return it
    def ensure_annotation(clade_el):
        ann = clade_el.find("px:annotation", NS)
        if ann is None:
            ann = ET.SubElement(clade_el, q("annotation"))
            ann.text = ""  # text container for GraPhlAn directives
        if ann.text is None:
            ann.text = ""
        return ann

    # Append a ring_color line (level 1) to the annotation text of a clade
    def add_ring_color_line(ann_el, color):
        # Keep existing text and ensure newline separation
        lines = (ann_el.text or "").splitlines()
        # If a ring_color line exists already, replace its color; else append a new line
        replaced = False
        new_lines = []
        for ln in lines:
            if "\tring_color\t" in ln:
                parts = ln.split("\t")
                if len(parts) >= 4:
                    parts[-1] = color
                    ln = "\t".join(parts)
                    replaced = True
            new_lines.append(ln)
        if not replaced:
            new_lines.append(f"ring_color\t1\t{color}")
        ann_el.text = "\n".join(l for l in new_lines if l.strip())

    # Patch loop
    patched_ok = 0
    for leaf, color in to_patch:
        cl = find_clade_by_name_fast(root, leaf)
        if cl is None:
            # leaf label mismatch? skip silently but report
            # (Your labels are already normalized consistently, so this should be rare.)
            # print(f"  ! Leaf not found in XML: {leaf}")
            continue
        ann = ensure_annotation(cl)
        before = (ann.text or "").count("ring_color")
        add_ring_color_line(ann, color.upper())
        after = (ann.text or "").count("ring_color")
        if after > before:
            patched_ok += 1

    # write out
    tree.write(XML_OUT, encoding="utf-8", xml_declaration=True)
    print(f"✅ Wrote patched XML: {XML_OUT} (added ring_color to {patched_ok} leaves)")

    # quick sanity counts
    txt_in  = XML_IN.read_text(errors="ignore")
    txt_out = XML_OUT.read_text(errors="ignore")
    print("ring_color count:", txt_in.count("ring_color"), "→", txt_out.count("ring_color"))
    leaf_cnt = sum(1 for _ in re.finditer(r"<name>[^<]+</name>", txt_out))
    print("leaf <name> tags in patched XML:", leaf_cnt)

    # render normal + big + svg
    def run(cmd):
        print("Running:", cmd)
        subprocess.run(shlex.split(cmd), check=True)

    run(f"graphlan.py {shlex.quote(str(XML_OUT))} {shlex.quote(str(PNG_OUT))} --dpi 300 --size 7.5 --pad 0.24 --avoid_reordering")
    run(f"graphlan.py {shlex.quote(str(XML_OUT))} {shlex.quote(str(BIG_PNG))} --dpi 600 --size 12 --pad 0.35 --avoid_reordering")
    run(f"graphlan.py {shlex.quote(str(XML_OUT))} {shlex.quote(str(SVG_OUT))} --format svg --size 12 --pad 0.35 --avoid_reordering")
    print("✅ Figures:", PNG_OUT, BIG_PNG, SVG_OUT)

donors = []
for m in donor_pat.finditer(txt):
    ann_block = m.group("ann")
    if "ring_color" in ann_block:
        donors.append((m.group("leaf"), ann_block))
        break

if donors:
    donor_leaf, donor_block = donors[0]
    print("Using donor leaf with ring block:", donor_leaf)
else:
    # Fallback: minimal block is enough because ring globals are already set
    print("No donor <annotation> with ring_color found — using minimal ring block fallback.")
    donor_block = "ring_color\t1\t#084594"  # placeholder color; will be recolored per-leaf
print("Any 'ring_color' text already in XML?", "ring_color" in txt)

# 4) Minimal understanding of the donor block:
#    It typically contains lines like: "ring_label\t1\t..." and "ring_color\t1\t#RRGGBB" etc.
#    We'll adjust only the color line(s) to the desired color for each target leaf,
#    and reuse the same ring level/structure that already renders on your GraPhlAn install.
#
#    Heuristic: replace every "#RRGGBB" in ring_color lines within donor_block with the target color.
def recolor_block(block_text: str, new_color: str) -> str:
    new_color = new_color.upper()
    changed = False
    out = []
    for ln in block_text.splitlines():
        if "\tring_color\t" in ln:
            parts = ln.split("\t")
            if len(parts) >= 4:
                parts[-1] = new_color
                ln = "\t".join(parts)
                changed = True
            else:
                # Replace any hex if present
                ln2 = re.sub(r"#([0-9A-Fa-f]{6})", new_color, ln)
                if ln2 != ln:
                    ln = ln2
                    changed = True
                else:
                    # No color present -> append it explicitly
                    ln = (ln.rstrip("\n") + "\t" + new_color)
                    changed = True
        out.append(ln)
    # If donor block had no ring_color line for some reason, add one
    if not any("\tring_color\t" in l for l in out):
        out.append(f"ring_color\t1\t{new_color}")
    return "\n".join(out)

# 5) Build a mapping: leaf -> target color from your *original* annotation text.
#    We want the intended ring_color for those skipped leaves.
ANNOT_TXT = Path("graphlan_annotations.txt")
if not ANNOT_TXT.exists():
    raise FileNotFoundError("Cannot find graphlan_annotations.txt to derive the intended ring colors.")

target_colors = {}  # leaf -> color
with ANNOT_TXT.open("r", encoding="utf-8") as fh:
    for raw in fh:
        raw = raw.strip()
        if not raw or raw.startswith("#"): 
            continue
        parts = raw.split("\t")
        if len(parts) == 4 and parts[1] == "ring_color":
            leaf, _, level, color = parts
            # keep first occurrence (level-1) or overwrite? We’ll keep the first.
            target_colors.setdefault(leaf, color)

missing_in_annot = [l for l in df["leaf"] if l not in target_colors]
if missing_in_annot:
    print(f"⚠️ {len(missing_in_annot)} skipped leaves have no ring_color in annot file (won't be patched). Example:", missing_in_annot[:5])


# 6) Inject a recolored donor block into each target leaf’s <clade>.
#    We’ll find the clade by <clade><name>LEAF</name> ... </clade>
#    Strategy:
#       - If the clade already has <annotation>...</annotation>, we will *append* a recolored ring_color line
#         to that block (but avoid duplicates).
#       - Else, we’ll create a fresh <annotation>...</annotation> using the recolored donor block.
#
#    We'll do this with regexes scoped to the specific clade to avoid breaking the XML.
def inject_ring_for_leaf(xml_text: str, leaf_name: str, color: str) -> str:
    # Match the specific clade section
    clade_re = re.compile(
        rf"(<clade>\s*<name>{re.escape(leaf_name)}</name>)(?P<body>.*?)(</clade>)",
        re.DOTALL
    )
    m = clade_re.search(xml_text)
    if not m:
        return xml_text  # leaf not found (shouldn't happen)
    head, body, tail = m.group(1), m.group("body"), m.group(3)

    if "<annotation>" in body:
        # Update or append a ring_color line inside existing annotation
        ann_re = re.compile(r"(<annotation>)(?P<ann>.*?)(</annotation>)", re.DOTALL)
        ma = ann_re.search(body)
        if not ma:
            # should not happen if "<annotation>" in body; but be safe
            new_ann = recolor_block(donor_block, color)
            new_body = body + f"\n<annotation>\n{new_ann}\n</annotation>\n"
        else:
            ann_txt = ma.group("ann")
            # If a ring_color for level 1 already exists, replace its color; else append a new line
            if "\tring_color\t" in ann_txt:
                # Replace color on all ring_color rows
                ann_txt2 = recolor_block(ann_txt, color)
            else:
                # Append one ring_color row (use donor format)
                donor_one = recolor_block(donor_block, color)
                # Append only lines that contain ring_color from donor to avoid duplicating headers
                donor_ring_lines = "\n".join([ln for ln in donor_one.splitlines() if "\tring_color\t" in ln])
                ann_txt2 = (ann_txt.rstrip() + "\n" + donor_ring_lines + "\n")
            new_body = body[:ma.start("ann")] + ann_txt2 + body[ma.end("ann"):]
    else:
        # No annotation block: create a new one from donor
        new_ann = recolor_block(donor_block, color)
        new_body = body + f"\n<annotation>\n{new_ann}\n</annotation>\n"

    # Replace the clade section
    start, end = m.span()
    return xml_text[:start] + head + new_body + tail + xml_text[end:]

# 7) Apply patches
to_patch = [(leaf, target_colors[leaf]) for leaf in df["leaf"] if leaf in target_colors]
print(f"Patching {len(to_patch)} leaves with ring_color…")
patched = txt
patched_count = 0
for leaf, color in to_patch:
    before = patched.count("ring_color")
    patched = inject_ring_for_leaf(patched, leaf, color)
    after = patched.count("ring_color")
    if after > before:
        patched_count += 1

print(f"Leaves that gained a ring_color entry: {patched_count} / {len(to_patch)}")

# 8) Write out and render
XML_OUT.write_text(patched, encoding="utf-8")
print("Wrote:", XML_OUT)

cmd = f"graphlan.py {shlex.quote(str(XML_OUT))} {shlex.quote(str(PNG_OUT))} --dpi 300 --size 7.5 --pad 0.24 --avoid_reordering"
print("Rendering:", cmd)
subprocess.run(shlex.split(cmd), check=True)
print("✅ Figure written:", PNG_OUT)

# Re-render the patched XML bigger + sharper + vector copy
BIG_PNG = Path("tree_circular_patched_big.png")
SVG_OUT = Path("tree_circular_patched.svg")

# Bigger canvas (inches) + higher DPI
cmd = f"graphlan.py {shlex.quote(str(XML_OUT))} {shlex.quote(str(BIG_PNG))} --dpi 600 --size 12 --pad 0.35 --avoid_reordering"
print("Rendering big PNG:", cmd)
subprocess.run(shlex.split(cmd), check=True)

# Also produce a vector SVG (great for zooming without blur)
cmd = f"graphlan.py {shlex.quote(str(XML_OUT))} {shlex.quote(str(SVG_OUT))} --format svg --size 12 --pad 0.35 --avoid_reordering"
print("Rendering SVG:", cmd)
subprocess.run(shlex.split(cmd), check=True)

print("✅ Wrote:", BIG_PNG, "and", SVG_OUT)

# Make ring taller and labels larger (applied to level 1 in-place)
BIGGER_XML = Path("taxonomy_tree_annotated.patched.bigger.xml")

# If you have _run_annot_with_lines from your diagnostic helpers, use it.
# Otherwise, here's a minimal inline helper:
def _run_annot_inline(in_xml: Path, lines, out_xml: Path) -> bool:
    import tempfile, shlex, subprocess
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as tmp:
        tmp.write("\n".join(lines) + "\n")
        tmp_path = Path(tmp.name)
    cmd = f"graphlan_annotate.py --annot {shlex.quote(str(tmp_path))} {shlex.quote(str(in_xml))} {shlex.quote(str(out_xml))}"
    res = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    ok = out_xml.exists() and out_xml.stat().st_size > 0
    if not ok:
        print("annotate failed:", res.returncode)
        print(res.stderr)
    return ok

# Tweak knobs here:
NEW_TITLE_SIZE = 20
NEW_RING_LABEL_SIZE = 13
NEW_RING_HEIGHT = 0.18    # was 0.12
NEW_RING_WIDTH  = 1.6     # was 1.0

bump_lines = [
    f"title_font_size\t{NEW_TITLE_SIZE}",
    f"ring_label_font_size\t1\t{NEW_RING_LABEL_SIZE}",
    f"ring_height\t1\t{NEW_RING_HEIGHT}",
    f"ring_width\t1\t{NEW_RING_WIDTH}",
]

print("Applying bigger ring/label/title settings…")
ok = _run_annot_inline(XML_OUT, bump_lines, BIGGER_XML)
print("ok?", ok)

# Render this “bigger” variant too:
BIGGER_PNG = Path("tree_circular_patched_bigger.png")
cmd = f"graphlan.py {shlex.quote(str(BIGGER_XML))} {shlex.quote(str(BIGGER_PNG))} --dpi 600 --size 12 --pad 0.35 --avoid_reordering"
print("Rendering bigger PNG:", cmd)
subprocess.run(shlex.split(cmd), check=True)
print("✅ Wrote:", BIGGER_PNG)
