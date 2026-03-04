import os, re, csv, shlex, subprocess, tempfile, shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import Counter
import pandas as pd
from io import StringIO
from Bio import Phylo

from pathlib import Path

proteome_dir = Path("../mycobacterium_files")
csv_files = list(proteome_dir.glob("*.csv"))
print(f"📦 Total .csv proteome files: {len(csv_files)}")

# ========= USER SETTINGS =========200
PROTEOME_CSV = Path("../mycobacterium_files/mycobacterium_proteomes.csv")
PROTEOME_DIR = Path("../mycobacterium_files")
TAXON_TERM   = "Mycobacterium"
BRANCH_LEN   = 1
RING_LEVEL   = 1
PNG_SIZE     = 7
PNG_DPI      = 300

# ========= HELPERS =========
CANONICAL_RANKS = ["domain","kingdom","phylum","class","order","family","genus","species"]
SAFE_CHARS_RE   = re.compile(r"[^A-Za-z0-9_.|\-]+")

def safe_label(name: str) -> str:
    if not isinstance(name, str) or not name.strip():
        return "unnamed"
    s = re.sub(r"[ /,():;]+", "_", name.strip())
    return re.sub(r"_+", "_", s).strip("_")

def split_lineage(lineage: str) -> List[str]:
    if not isinstance(lineage, str) or not lineage.strip():
        return []
    sep = "," if "," in lineage else ";"
    return [p.strip() for p in lineage.split(sep) if p.strip()]

def clean_taxon_name(tag: str) -> str:
    if not isinstance(tag, str):
        return ""
    m = re.match(r"^(.*?)(?:\s*\([^)]*\))?$", tag.strip())
    name = m.group(1) if m else tag
    return re.sub(r"\s+", " ", name.replace('"', "").replace("'", "")).strip()

def lineage_to_nodes(lineage: str) -> List[str]:
    return [safe_label(clean_taxon_name(p)) for p in split_lineage(lineage) if p.strip()]

def extract_species_from_lineage(lineage: str) -> str:
    """Extract the final species name (rank = species) from a lineage string."""
    parts = split_lineage(lineage)
    for p in reversed(parts):
        if "(species)" in p.lower():
            return safe_label(clean_taxon_name(p.replace("(species)", "")))
    # fallback: use last token
    return safe_label(clean_taxon_name(parts[-1] if parts else "unnamed"))

def trim_lineage_to_term(lineage: str, term: str) -> Optional[str]:
    if not lineage:
        return None
    parts = split_lineage(lineage)
    clean_parts = [clean_taxon_name(p).lower() for p in parts]
    term_l = term.lower()
    for i, p in enumerate(clean_parts):
        if term_l in p:
            return "; ".join(parts[i:])
    return None

def write_graphlan_annotations(
    leaf2pct: Dict[str, float],
    out_path: Path,
    ring_level: int = 1
):
    lines = [
        "title\tSAMs across species",
        "branch_thickness\t0.7",
        "branch_color\t#666666",
        "total_plotted_degrees\t340",
        f"ring_label\t{ring_level}\t% SAMs",
        f"ring_label_font_size\t{ring_level}\t11",
        f"ring_internal_separator_thickness\t{ring_level}\t0.0",
        f"ring_external_separator_thickness\t{ring_level}\t0.0",
        f"ring_separator_color\t{ring_level}\t#DDDDDD",
        f"ring_width\t{ring_level}\t1.0",
        f"ring_height\t{ring_level}\t0.25",
    ]

    for leaf, v in leaf2pct.items():
        color = value_to_hex(v, 0.0, 1.0)
        lines.append(f"{leaf}\tclade_marker_color\t{color}")
        lines.append(f"{leaf}\tclade_marker_size\t36")
        lines.append(f"{leaf}\tclade_marker_shape\to")  # circle
        lines.append(f"{leaf}\tannotation\t{leaf.replace('unclassified_', 'u.')}")
        lines.append(f"{leaf}\tannotation_font_size\t12")
        lines.append(f"{leaf}\tannotation_background_color\t#FFFFFF")
        lines.append(f"{leaf}\tring_color\t{ring_level}\t{color}")

    print("🧪 Validating annotation file contents...")
    seen_classes = 0
    malformed = 0
    unknown_targets = 0
    clade_names = set(leaf2pct.keys())

    for line in lines:
        parts = line.split("\t")
        if len(parts) == 0:
            continue

        if len(parts) not in {2, 3, 4}:
            print(f"❌ Malformed line (wrong number of fields): {line}")
            malformed += 1

        if len(parts) >= 2 and parts[1].strip().lower() == "classes":
            print(f"❌ ILLEGAL keyword 'classes' found: {line}")
            seen_classes += 1

        # Catch lines targeting nonexistent clades
        if len(parts) >= 1 and parts[1] not in {
            "title", "branch_thickness", "branch_color", "total_plotted_degrees"
        } and not parts[0].startswith("*") and parts[0] not in clade_names:
            print(f"⚠️ Targeted clade not found: {parts[0]}")
            unknown_targets += 1

    print(f"✅ Annotation validation complete.")
    print(f"    → Illegal 'classes' lines: {seen_classes}")
    print(f"    → Malformed lines:         {malformed}")
    print(f"    → Unknown clade targets:   {unknown_targets}")

    # Save annotation file
    out_path.write_text("\n".join(lines) + "\n")
    print(f"📄 Annotation file written to: {out_path}")

def value_to_hex(v: float, vmin=0.0, vmax=1.0) -> str:
    C0, C1 = "#e9eef3", "#084594"
    if v is None:
        return "#BDBDBD"
    v = min(max(v, vmin), vmax)
    t = 0.0 if vmax == vmin else (v - vmin)/(vmax - vmin)
    r0,g0,b0 = (int(C0[1:3],16), int(C0[3:5],16), int(C0[5:7],16))
    r1,g1,b1 = (int(C1[1:3],16), int(C1[3:5],16), int(C1[5:7],16))
    r = int(round(r0 + (r1 - r0)*t))
    g = int(round(g0 + (g1 - g0)*t))
    b = int(round(b0 + (b1 - b0)*t))
    return f"#{r:02X}{g:02X}{b:02X}"

# ========= TREE STRUCTURE =========
class TreeNode:
    def __init__(self, name=None):
        self.name = name
        self.children: Dict[str, "TreeNode"] = {}

    def add_path(self, path: List[str], leaf_label: Optional[str]=None):
        node = self
        for n in path:
            node = node.children.setdefault(n, TreeNode(n))
        if leaf_label:
            node.children.setdefault(leaf_label, TreeNode(leaf_label))

def to_newick(node: TreeNode, is_root=False, branch_len=1) -> str:
    if not node.children:
        return f"{node.name}:{branch_len}"
    subs = ",".join(to_newick(c, False, branch_len) for c in node.children.values())
    label = node.name if node.name else ""
    return f"({subs}){label}" if is_root else f"({subs}){label}:{branch_len}"

# ========= CORE STEPS =========
from pathlib import Path
import pandas as pd

def build_sam_percentages_and_lineages(metadata_df: pd.DataFrame, proteome_dir: Path):
    sam_data = []
    lineage_map = {}

    # Clean organism names from metadata for matching
    metadata_df["Organism_clean"] = metadata_df["Organism"].str.replace("Mycobacterium", "").str.strip().str.replace(r"\(.*\)", "", regex=True).str.replace(" ", "_").str.lower()

    for proteome_file in proteome_dir.glob("*.csv"):
        stem_name = proteome_file.stem.lower()

        # Match file name to metadata
        matched_row = metadata_df[metadata_df["Organism_clean"].apply(lambda x: x in stem_name)]
        if matched_row.empty:
            print(f"❌ No match found in metadata for: {proteome_file.name}")
            continue

        row = matched_row.iloc[0]
        lineage_raw = row["Taxonomic lineage"]
        lineage_parts = [x.strip().replace(" ", "_") for x in lineage_raw.split(",") if x.strip()]
        leaf_label = proteome_file.stem.replace("_", "_")  # use full filename without extension
        lineage_parts.append(leaf_label)

        try:
            df = pd.read_csv(proteome_file)
        except Exception as e:
            print(f"❌ Could not read {proteome_file.name}: {e}")
            continue

        if "InterPro" not in df.columns:
            print(f"⚠️ Missing 'InterPro' column in {proteome_file.name}")
            continue

        total_proteins = len(df)
        if "InterPro" not in df.columns:
            print(f"⚠️ Missing 'InterPro' column in {proteome_file.name}")
            continue

        print(f"🔍 InterPro dtype in {proteome_file.name}: {df['InterPro'].dtype}")

        interpro_series = df["InterPro"].astype(str)
        sam_hits = interpro_series.str.contains("IPR029063", na=False).sum()
        percent_sam = (sam_hits / total_proteins) * 100 if total_proteins else 0
        sam_data.append((leaf_label, percent_sam))
        lineage_map[leaf_label] = lineage_parts

        print(f"✅ {leaf_label}: {sam_hits}/{total_proteins} SAMs ({percent_sam:.2f}%)")

    return sam_data, lineage_map

def annotate_and_render(xml_in: Path, annot_txt: Path, xml_out: Path, png_out: Path):
    ann_bin = shutil.which("graphlan_annotate.py")
    if ann_bin is None:
        raise RuntimeError("graphlan_annotate.py not found in PATH!")
    png_bin = shutil.which("graphlan.py")
    if png_bin is None:
        raise RuntimeError("graphlan.py not found in PATH!")

    cmd_ann = f"{ann_bin} --annot {annot_txt} {xml_in} {xml_out}"
    print(f"🔧 Running annotation command: {cmd_ann}")
    if not xml_in.exists():
        raise FileNotFoundError(f"❌ Input XML file not found: {xml_in}")
    else:
        print(f"📁 Input XML file exists: {xml_in}")
    ann = subprocess.run(shlex.split(cmd_ann), capture_output=True, text=True)
    if ann.returncode != 0:
        print("❌ Annotation command failed!")
        print("stdout:", ann.stdout)
        print("stderr:", ann.stderr)
        raise RuntimeError("graphlan_annotate.py failed.")
    else:
        print("✅ Annotation command completed successfully.")
    if ann.returncode != 0:
        raise RuntimeError("graphlan_annotate.py failed.")
    print("stdout:", ann.stdout)
    print("stderr:", ann.stderr)

    # New debug check
    if not xml_out.exists():
        print(f"❌ Annotated file was not created: {xml_out}")
        print("📦 Contents of working directory:", list(Path(".").glob("*")))
        raise RuntimeError("Annotated XML file is missing after annotation step!")
    else:
        print(f"✅ Annotated file created: {xml_out}")
    cmd_png = f"{png_bin} {xml_out} {png_out} --dpi {PNG_DPI} --size {PNG_SIZE} --pad 0.24 --avoid_reordering"
    print(f"🎨 Running render command: {cmd_png}")
    png = subprocess.run(shlex.split(cmd_png), capture_output=True, text=True)
    if png.returncode != 0:
        raise RuntimeError("graphlan.py failed.")
    print("stdout:", png.stdout)
    print("stderr:", png.stderr)

def build_tree(pairs: List[Tuple[str,str]]) -> Tuple[TreeNode, str]:
    root, dup = TreeNode(None), Counter()
    for leaf_label, lineage in pairs:
        path = lineage_to_nodes(lineage)
        if not path:
            print(f"⚠️ No valid node path for {leaf_label} → Skipping")
            continue

        dup[leaf_label] += 1
        # Make sure the leaf label is unique and doesn't conflict with internal names
        clean_leaf = leaf_label
        if dup[leaf_label] > 1 or clean_leaf in path:
            clean_leaf = f"{leaf_label}__{dup[leaf_label]}"
        print(f"🌱 Added leaf {clean_leaf} under path: {' > '.join(path)}")
        node = root
        for n in path:
            node = node.children.setdefault(n, TreeNode(n))
        node.children.setdefault(clean_leaf, TreeNode(clean_leaf))

    newick = to_newick(root, is_root=True, branch_len=BRANCH_LEN) + ";"
    return root, newick
# ========= MAIN =========
def main():
    out_prefix = TAXON_TERM.replace(" ", "_")
    out_xml     = Path(f"{out_prefix}_tree.xml")
    out_xml_ann = Path(f"{out_prefix}_annotated.xml")
    out_png     = Path(f"{out_prefix}_tree.png")
    out_annot   = Path(f"{out_prefix}_annotations.txt")

    # Build lineage and %SAMs from the proteome CSVs
    df_meta = pd.read_csv(PROTEOME_CSV, dtype=str)
    leaf2pct, pairs = build_sam_percentages_and_lineages(df_meta, PROTEOME_DIR)
    if not pairs:
        raise SystemExit("❌ No valid lineage/species pairs found.")

    df_trim = pd.DataFrame(pairs, columns=["leaf_label", "lineage"]).drop_duplicates(subset=["leaf_label"])
    pairs = list(zip(df_trim["leaf_label"], df_trim["lineage"]))
    print(f"🧬 Number of label-lineage pairs before deduplication: {len(pairs)}")
    df_trim = df_trim.drop_duplicates(subset=["leaf_label"])
    print(f"🧬 Number of unique leaf labels: {df_trim['leaf_label'].nunique()}")
    print(f"🧪 Tree build input sample (first 5 pairs):")
    for leaf, lineage in pairs[:5]:
        print(f" - {leaf}: {lineage}")
    empty_path_count = 0
    for leaf, lineage in pairs:
        path = lineage_to_nodes(lineage)
        if not path:
            print(f"⚠️ Lineage produced empty node path: {lineage}")
            empty_path_count += 1
    print(f"📉 {empty_path_count}/{len(pairs)} pairs produced empty paths")
    leaf_counts = Counter([leaf for leaf, _ in pairs])
    dupes = {k: v for k, v in leaf_counts.items() if v > 1}
    if dupes:
        print("⚠️ Duplicate leaf labels detected:")
        for name, count in dupes.items():
            print(f" - {name}: {count} occurrences")
    tree_root, newick = build_tree(pairs)
    tree = Phylo.read(StringIO(newick), "newick")
    internal_nodes = [cl for cl in tree.get_nonterminals() if cl.name]
    Phylo.write(tree, out_xml, "phyloxml")
    leaves = [cl.name for cl in tree.get_terminals() if cl.name]
    if not leaves:
        print("❌ No leaves found in tree! Tree construction likely failed.")
        print("📉 Check if lineage parsing is failing or leaf labels are incorrect.")
        return
    internal_names = [cl.name for cl in internal_nodes if cl.name not in leaves]
    print(f"🌳 Tree written to {out_xml}")
    print(f"🌿 Tree leaves count: {len(leaves)}")
    print(f"🧩 Internal nodes count: {len(internal_names)}")
    # Collect RNA-binding % using consistent species names
    # Match only tree leaves
    leaf2pct = {leaf: leaf2pct.get(leaf, 0.0) for leaf in leaves}
    write_graphlan_annotations(leaf2pct, out_annot, ring_level=RING_LEVEL)
    print(f"📄 Wrote annotation file: {out_annot}")
    lines = out_annot.read_text().strip().splitlines()
    print(f"📝 Annotation file has {len(lines)} lines.")
    has_clades = any("\tclade_label\t" in line for line in lines)
    if not has_clades:
        print("⚠️ No clade labels found in the annotation file.")
    try:
        tree_test = Phylo.read(out_xml, "phyloxml")
        print("✅ Successfully re-read written PhyloXML file.")
    except Exception as e:
        print("❌ Failed to parse the written PhyloXML file:")
        print(e)
    annotate_and_render(out_xml, out_annot, out_xml_ann, out_png)
    if not out_xml_ann.exists():
        print(f"❌ Annotated XML file not created: {out_xml_ann}")
    else:
        print(f"✅ Annotated XML created: {out_xml_ann}")
    if not out_png.exists():
        print(f"❌ PNG file not created: {out_png}")
    else:
        print(f"✅ PNG output created at: {out_png}")
    print("🧬 Annotating the following leaves:")
    for leaf in leaf2pct:
        print(f"  - {leaf}")
    tree_leaves = [cl.name for cl in tree.get_terminals()]
    ann_leaves = set(leaf2pct.keys())
    missing_in_tree = ann_leaves - set(tree_leaves)

    if missing_in_tree:
        print("⚠️ Warning: Some annotated leaves do not appear in the tree:")
        for name in missing_in_tree:
            print(f"  - {name}")
    print(f"🧬 Tree has {len(tree.get_terminals())} leaves.")
    if len(tree.get_terminals()) < 2:
        print("⚠️ Tree is too small — this might prevent GraPhlAn from rendering.")
    print("🧬 Tree leaves:")
    for leaf in tree.get_terminals():
        print(f" - {leaf.name}")

    print("🧬 Annotated leaves:")
    for leaf in leaf2pct:
        print(f" - {leaf}")
    missing_in_tree = set(leaf2pct) - {cl.name for cl in tree.get_terminals()}
    if missing_in_tree:
        print("⚠️ WARNING: Some annotated leaves are not found in the tree:")
        for leaf in missing_in_tree:
            print(f"  - {leaf}")

    print(f"✅ Final tree has {len(tree.get_terminals())} leaf nodes.")
    print(f"✅ Final annotated leaf2pct contains {len(leaf2pct)} entries.")


if __name__ == "__main__":
    main()