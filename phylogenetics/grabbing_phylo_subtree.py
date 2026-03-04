# import os, re, csv, shlex, subprocess, tempfile, shutil
# from pathlib import Path
# from typing import Dict, List, Tuple, Optional
# from collections import Counter
# import pandas as pd
# from io import StringIO
# from Bio import Phylo

# # ========= USER SETTINGS =========200
# PRED_DIR     = Path("../prediction_files")
# SCORE1_PATH  = Path("../Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03.tsv.gz")
# SCORE4_PATH  = Path("../Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03 (1).tsv.gz")
# TAXON_TERM   = "Bacteria"
# BRANCH_LEN   = 1
# RING_LEVEL   = 1
# PNG_SIZE     = 7
# PNG_DPI      = 300

# # ========= HELPERS =========
# CANONICAL_RANKS = ["domain","kingdom","phylum","class","order","family","genus","species"]
# SAFE_CHARS_RE   = re.compile(r"[^A-Za-z0-9_.|\-]+")

# def safe_label(name: str) -> str:
#     if not isinstance(name, str) or not name.strip():
#         return "unnamed"
#     s = re.sub(r"[ /,():;]+", "_", name.strip())
#     return re.sub(r"_+", "_", s).strip("_")

# def split_lineage(lineage: str) -> List[str]:
#     if not isinstance(lineage, str) or not lineage.strip():
#         return []
#     sep = "," if lineage.count(",") > lineage.count(";") else ";"
#     return [p.strip() for p in lineage.split(sep) if p.strip()]

# def clean_taxon_name(tag: str) -> str:
#     if not isinstance(tag, str):
#         return ""
#     m = re.match(r"^(.*?)(?:\s*\([^)]*\))?$", tag.strip())
#     name = m.group(1) if m else tag
#     return re.sub(r"\s+", " ", name.replace('"', "").replace("'", "")).strip()

# def lineage_to_nodes(lineage: str) -> List[str]:
#     nodes = []
#     for raw in split_lineage(lineage):
#         m = re.match(r"^(.*?)(?:\s*\(([^)]+)\))?$", raw)
#         if not m:
#             continue
#         name, rank = clean_taxon_name(m.group(1)), (m.group(2) or "").lower()
#         if name and rank in CANONICAL_RANKS and rank != "species":
#             nodes.append(name)
#     return nodes

# def extract_species_from_lineage(lineage: str) -> str:
#     """Extract the final species name (rank = species) from a lineage string."""
#     parts = split_lineage(lineage)
#     for p in reversed(parts):
#         if "(species)" in p.lower():
#             return safe_label(clean_taxon_name(p.replace("(species)", "")))
#     # fallback: use last token
#     return safe_label(clean_taxon_name(parts[-1] if parts else "unnamed"))

# def trim_lineage_to_term(lineage: str, term: str) -> Optional[str]:
#     if not lineage:
#         return None
#     parts = split_lineage(lineage)
#     clean_parts = [clean_taxon_name(p).lower() for p in parts]
#     term_l = term.lower()
#     for i, p in enumerate(clean_parts):
#         if term_l in p:
#             return "; ".join(parts[i:])
#     return None

# def value_to_hex(v: float, vmin=0.0, vmax=1.0) -> str:
#     C0, C1 = "#e9eef3", "#084594"
#     if v is None:
#         return "#BDBDBD"
#     v = min(max(v, vmin), vmax)
#     t = 0.0 if vmax == vmin else (v - vmin)/(vmax - vmin)
#     r0,g0,b0 = (int(C0[1:3],16), int(C0[3:5],16), int(C0[5:7],16))
#     r1,g1,b1 = (int(C1[1:3],16), int(C1[3:5],16), int(C1[5:7],16))
#     r = int(round(r0 + (r1 - r0)*t))
#     g = int(round(g0 + (g1 - g0)*t))
#     b = int(round(b0 + (b1 - b0)*t))
#     return f"#{r:02X}{g:02X}{b:02X}"

# # ========= TREE STRUCTURE =========
# class TreeNode:
#     def __init__(self, name=None):
#         self.name = name
#         self.children: Dict[str, "TreeNode"] = {}

#     def add_path(self, path: List[str], leaf_label: Optional[str]=None):
#         node = self
#         for n in path:
#             node = node.children.setdefault(n, TreeNode(n))
#         if leaf_label:
#             node.children.setdefault(leaf_label, TreeNode(leaf_label))

# def to_newick(node: TreeNode, is_root=False, branch_len=1) -> str:
#     if not node.children:
#         return f"{node.name}:{branch_len}"
#     subs = ",".join(to_newick(c, False, branch_len) for c in node.children.values())
#     label = node.name if node.name else ""
#     return f"({subs}){label}" if is_root else f"({subs}){label}:{branch_len}"

# # ========= CORE STEPS =========
# def build_species_lineages(pred_dir: Path, score1_path: Path, score4_path: Path) -> pd.DataFrame:
#     score1 = pd.read_csv(score1_path, sep="\t", compression="gzip", dtype=str)
#     score4 = pd.read_csv(score4_path, sep="\t", compression="gzip", dtype=str)
#     entry2lin = {**dict(zip(score4["Entry"], score4["Taxonomic lineage"])),
#                  **dict(zip(score1["Entry"], score1["Taxonomic lineage"]))}  # score1 overrides

#     species2lineages: Dict[str, List[str]] = {}

#     for p in pred_dir.glob("*.csv"):
#         try:
#             df = pd.read_csv(p, dtype=str, encoding_errors="ignore", on_bad_lines="skip")
#         except Exception:
#             continue
#         if not {"Organism", "Entry"}.issubset(df.columns):
#             continue
#         for _, row in df.iterrows():
#             org = safe_label(row["Organism"])
#             entry = row["Entry"]
#             if pd.notna(entry) and entry in entry2lin:
#                 species2lineages.setdefault(org, []).append(entry2lin[entry])

#     rows = []
#     for species, lin_list in species2lineages.items():
#         # Pick the most detailed lineage
#         lin_list.sort(key=lambda x: len(split_lineage(x)), reverse=True)
#         rows.append({"leaf_label": species, "lineage": lin_list[0]})
#     return pd.DataFrame(rows), entry2lin

# def collect_rna_values(pred_dir: Path, entry2lin: Dict[str, str], taxon_term: str, thresh=0.5) -> Dict[str, float]:
#     """Return {species → RNA-binding %}, using Entry to get consistent species label via lineage."""
#     from collections import defaultdict

#     species_vals: Dict[str, List[float]] = defaultdict(list)

#     for f in pred_dir.glob("*.csv"):
#         try:
#             df = pd.read_csv(f, dtype=str, encoding_errors="ignore", on_bad_lines="skip")
#         except Exception:
#             continue
#         if not {"Entry", "RNA", "Organism"}.issubset(df.columns):
#             continue
#         df = df.dropna(subset=["Entry", "RNA"])
#         df["RNA"] = pd.to_numeric(df["RNA"], errors="coerce")
#         for _, row in df.iterrows():
#             entry = row["Entry"]
#             rna_val = row["RNA"]
#             if pd.isna(rna_val) or entry not in entry2lin:
#                 continue
#             lineage = entry2lin[entry]
#             trimmed = trim_lineage_to_term(lineage, taxon_term)
#             if not trimmed:
#                 continue
#             species = extract_species_from_lineage(trimmed)
#             val = 1.0 if rna_val >= thresh else 0.0
#             species_vals[species].append(val)

#     return {k: float(sum(v)/len(v)) for k, v in species_vals.items() if v}

# def build_tree(pairs: List[Tuple[str,str]]) -> Tuple[TreeNode, str]:
#     root, dup = TreeNode(None), Counter()
#     for leaf_label, lineage in pairs:
#         path = lineage_to_nodes(lineage)
#         if not path:
#             continue

#         dup[leaf_label] += 1
#         # Make sure the leaf label is unique and doesn't conflict with internal names
#         clean_leaf = leaf_label
#         if dup[leaf_label] > 1 or clean_leaf in path:
#             clean_leaf = f"{leaf_label}__{dup[leaf_label]}"

#         node = root
#         for n in path:
#             node = node.children.setdefault(n, TreeNode(n))
#         node.children.setdefault(clean_leaf, TreeNode(clean_leaf))

#     newick = to_newick(root, is_root=True, branch_len=BRANCH_LEN) + ";"
#     return root, newick

# # ========= NEW: INTERACTIVE CYTOSCAPE-STYLE TREE USING PYVIS =========
# from pyvis.network import Network
# from Bio import Phylo

# def build_pyvis_tree(newick_str, leaf2pct, out_html="rna_tree.html"):
#     """
#     Build a Cytoscape-style interactive circular tree using PyVis (Vis.js).
#     - newick_str: string, Newick format tree
#     - leaf2pct: dict {leaf -> RNA-binding %}
#     - out_html: output filename
#     """

#     # 1. Parse the tree using BioPython
#     tree = Phylo.read(StringIO(newick_str), "newick")

#     # 2. Initialize PyVis network
#     net = Network(
#         height="900px",
#         width="100%",
#         bgcolor="#ffffff",
#         font_color="black",
#         directed=False
#     )

#     # Enable physics for circular layout
#     net.barnes_hut()
#     net.set_options("""
#     {
#       "layout": {
#         "improvedLayout": true
#       },
#       "physics": {
#         "enabled": true,
#         "barnesHut": {
#           "gravitationalConstant": -8000,
#           "centralGravity": 0.3,
#           "springLength": 100,
#           "springConstant": 0.04,
#           "avoidOverlap": 1
#         }
#       },
#       "nodes": {
#         "shape": "dot",
#         "scaling": { "min": 5, "max": 30 }
#       },
#       "edges": {
#         "smooth": true
#       }
#     }
#     """)

#     # 3. Traverse tree and add nodes/edges
#     def add_clade(clade, parent=None):
#         node_id = id(clade)

#         label = clade.name if clade.name else ""
#         size = 12

#         # Get color based on RNA-binding %
#         if label in leaf2pct:
#             pct = leaf2pct[label]
#             color = value_to_hex(pct)     # your gradient function
#             size = 20 + pct * 20          # bigger for high RNA-binding
#             title = f"{label}<br>RNA-binding: {pct:.2f}"
#         else:
#             color = "#D3D3D3"
#             title = label

#         # Add node to the graph
#         net.add_node(
#             node_id,
#             label=label,
#             title=title,
#             color=color,
#             size=size
#         )

#         if parent is not None:
#             net.add_edge(parent, node_id)

#         for child in clade.clades:
#             add_clade(child, parent=node_id)

#     # Start recursion
#     add_clade(tree.root)

#     # 4. Save interactive HTML
#     net.write_html(out_html)
#     print(f"🌐 Saved interactive Cytoscape-style tree to: {out_html}")

# def write_itol_gradient(leaf2pct, outfile="rna_binding_dataset.txt",
#                         dataset_label="RNA_binding",
#                         field_label="RNA %",
#                         legend_title="% RNA-binding",
#                         color="#084594"):
#     """
#     Write an iTOL DATASET_GRADIENT file for percent RNA-binding per leaf.

#     Parameters:
#     - leaf2pct: dict of {leaf_name: percent_value (0-1 or 0-100)}
#     - outfile: output filename (string)
#     - dataset_label: label shown in iTOL sidebar
#     - field_label: label shown in the gradient legend
#     - legend_title: title of the legend
#     - color: base color for the gradient (hex string)
#     """

#     with open(outfile, "w") as f:
#         f.write("DATASET_GRADIENT\n")
#         f.write("SEPARATOR TAB\n")
#         f.write(f"DATASET_LABEL\t{dataset_label}\n")
#         f.write(f"FIELD_LABEL\t{field_label}\n")
#         f.write(f"COLOR\t{color}\n")
#         f.write(f"LEGEND_TITLE\t{legend_title}\n")
#         f.write("DATA\n")
#         for leaf, pct in leaf2pct.items():
#             f.write(f"{leaf}\t{pct:.3f}\n")

# # ========= MAIN =========
# def main():
#     out_prefix = TAXON_TERM.replace(" ", "_")
#     out_xml     = Path(f"{out_prefix}_tree.xml")
#     out_xml_ann = Path(f"{out_prefix}_annotated.xml")
#     out_png     = Path(f"{out_prefix}_tree.png")
#     out_annot   = Path(f"{out_prefix}_annotations.txt")

#     df_used, entry2lin = build_species_lineages(PRED_DIR, SCORE1_PATH, SCORE4_PATH)
#     if df_used.empty:
#         raise SystemExit("No lineages found.")

#     keep_rows = []
#     for _, row in df_used.iterrows():
#         trimmed = trim_lineage_to_term(row["lineage"], TAXON_TERM)
#         if trimmed:
#             species = extract_species_from_lineage(row["lineage"])
#             keep_rows.append({"leaf_label": species, "lineage": trimmed})
#     df_trim = pd.DataFrame(keep_rows)
#     if df_trim.empty:
#         raise SystemExit(f"No lineage contains '{TAXON_TERM}'.")

#     pairs = list(zip(df_trim["leaf_label"], df_trim["lineage"]))
#     df_trim = df_trim.drop_duplicates(subset=["leaf_label"])
#     tree_root, newick = build_tree(pairs)
#     tree = Phylo.read(StringIO(newick), "newick")
#     internal_nodes = [cl for cl in tree.get_nonterminals() if cl.name]
#     Phylo.write(tree, out_xml, "phyloxml")
#     leaves = [cl.name for cl in tree.get_terminals() if cl.name]
#     internal_names = [cl.name for cl in internal_nodes if cl.name not in leaves]
#     # Collect RNA-binding % using consistent species names
#     leaf2pct = collect_rna_values(PRED_DIR, entry2lin, TAXON_TERM, thresh=0.5)
#     # Match only tree leaves
#     leaf2pct = {leaf: leaf2pct.get(leaf, 0.0) for leaf in leaves}
#     # build_pyvis_tree(newick, leaf2pct, out_html=f"{out_prefix}_rna_tree.html")
#     write_itol_gradient(leaf2pct)

# if __name__ == "__main__":
#     main()



# import os, re, csv, shlex, subprocess, tempfile, shutil
# from typing import Dict, List, Tuple, Optional
# from collections import Counter
# import pandas as pd
# from io import StringIO
# from Bio import Phylo
# from pathlib import Path

# # ========== USER INPUT ==========
# PRED_DIR     = Path("../mycobacterium_files")
# SCORE1_PATH  = Path("../Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03.tsv.gz")
# SCORE4_PATH  = Path("../Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03 (1).tsv.gz")
# BRANCH_LEN   = 1


# # ========== HELPERS ==========
# def safe_label(name: str) -> str:
#     if not isinstance(name, str) or not name.strip():
#         return "unnamed"
#     s = re.sub(r"[ /,():;]+", "_", name.strip())
#     return re.sub(r"_+", "_", s).strip("_")

# def split_lineage(lineage: str) -> List[str]:
#     if not isinstance(lineage, str) or not lineage.strip():
#         return []
#     sep = "," if lineage.count(",") > lineage.count(";") else ";"
#     return [p.strip() for p in lineage.split(sep) if p.strip()]

# def clean_taxon_name(tag: str) -> str:
#     if not isinstance(tag, str):
#         return ""
#     m = re.match(r"^(.*?)(?:\s*\([^)]*\))?$", tag.strip())
#     name = m.group(1) if m else tag
#     return re.sub(r"\s+", " ", name.replace('"', "").replace("'", "")).strip()

# def extract_species_from_lineage(lineage: str) -> str:
#     parts = split_lineage(lineage)
#     return safe_label(clean_taxon_name(parts[-1] if parts else "unnamed"))

# def lineage_to_nodes_all(lineage: str) -> List[str]:
#     return [safe_label(clean_taxon_name(p)) for p in split_lineage(lineage)]

# # ========== TREE STRUCTURE ==========
# class TreeNode:
#     def __init__(self, name=None):
#         self.name = name
#         self.children: Dict[str, "TreeNode"] = {}

#     def add_path(self, path: List[str], leaf_label: Optional[str]=None):
#         node = self
#         for n in path:
#             node = node.children.setdefault(n, TreeNode(n))
#         if leaf_label:
#             node.children.setdefault(leaf_label, TreeNode(leaf_label))

# def to_newick(node: TreeNode, is_root=False, branch_len=1) -> str:
#     if not node.children:
#         return f"{node.name}:{branch_len}"
#     subs = ",".join(to_newick(c, False, branch_len) for c in node.children.values())
#     label = node.name if node.name else ""
#     return f"({subs}){label}" if is_root else f"({subs}){label}:{branch_len}"

# # ========== MAIN WORKFLOW ==========
# def build_lineage_tree(pred_dir: Path, score1_path: Path, score4_path: Path, branch_len=1) -> Tuple[str, pd.DataFrame]:
#     score1 = pd.read_csv(score1_path, sep="\t", compression="gzip", dtype=str)
#     score4 = pd.read_csv(score4_path, sep="\t", compression="gzip", dtype=str)
#     entry2lin = {**dict(zip(score4["Entry"], score4["Taxonomic lineage"])),
#                  **dict(zip(score1["Entry"], score1["Taxonomic lineage"]))}  # score1 overrides

#     species2lineages: Dict[str, List[str]] = {}

#     for p in pred_dir.glob("*.csv"):
#         try:
#             df = pd.read_csv(p, dtype=str, encoding_errors="ignore", on_bad_lines="skip")
#         except Exception:
#             continue
#         if not {"Organism", "Entry"}.issubset(df.columns):
#             continue
#         for _, row in df.iterrows():
#             org = safe_label(row["Organism"])
#             entry = row["Entry"]
#             if pd.notna(entry) and entry in entry2lin:
#                 species2lineages.setdefault(org, []).append(entry2lin[entry])

#     rows = []
#     for species, lin_list in species2lineages.items():
#         lin_list.sort(key=lambda x: len(split_lineage(x)), reverse=True)
#         best_lineage = lin_list[0]
#         rows.append({"leaf_label": species, "lineage": best_lineage})

#     df_lineage = pd.DataFrame(rows)

#     # Build tree
#     root = TreeNode()
#     for _, row in df_lineage.iterrows():
#         lineage_nodes = lineage_to_nodes_all(row["lineage"])
#         root.add_path(lineage_nodes, leaf_label=row["leaf_label"])

#     newick = to_newick(root, is_root=True, branch_len=branch_len) + ";"
#     return newick, df_lineage

# tree_newick, df = build_lineage_tree(PRED_DIR, SCORE1_PATH, SCORE4_PATH)
# print(tree_newick[:500])  # preview Newick string
# # Save outputs
# out_dir = Path("./output_tree")  # or change this to any preferred path
# out_dir.mkdir(exist_ok=True)

# # Write Newick tree
# newick_path = out_dir / "myco_tree.nwk"
# with open(newick_path, "w") as f:
#     f.write(tree_newick)

# # Write PhyloXML
# tree = Phylo.read(StringIO(tree_newick), "newick")
# phyloxml_path = out_dir / "myco_tree.xml"
# Phylo.write(tree, phyloxml_path, "phyloxml")

# # Save lineage mapping table
# csv_path = out_dir / "myco_lineages.csv"
# df.to_csv(csv_path, index=False)

# print(f"✅ Tree written to: {newick_path}")
# print(f"✅ PhyloXML written to: {phyloxml_path}")
# print(f"✅ Lineage table written to: {csv_path}")


import pandas as pd
import re
from pathlib import Path
from typing import List, Dict, Optional
from collections import defaultdict

# === HELPERS ===
def safe_label(name: str) -> str:
    if not isinstance(name, str) or not name.strip():
        return "unnamed"
    s = re.sub(r"[ /,():;]+", "_", name.strip())
    return re.sub(r"_+", "_", s).strip("_")

def split_lineage(lineage: str) -> List[str]:
    if not isinstance(lineage, str) or not lineage.strip():
        return []
    sep = "," if lineage.count(",") > lineage.count(";") else ";"
    return [p.strip() for p in lineage.split(sep) if p.strip()]

def clean_taxon_name(tag: str) -> str:
    if not isinstance(tag, str):
        return ""
    m = re.match(r"^(.*?)(?:\s*\([^)]*\))?$", tag.strip())
    name = m.group(1) if m else tag
    return re.sub(r"\s+", " ", name.replace('"', "").replace("'", "")).strip()

def lineage_to_nodes_all(lineage: str) -> List[str]:
    return [safe_label(clean_taxon_name(p)) for p in split_lineage(lineage)]

# === TREE STRUCTURE ===
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

# === ITOL COLORSTRIP WRITER ===
def write_itol_colorstrip(data: Dict[str, float], out_file: Path, label: str, hex_color: str):
    with open(out_file, "w") as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"DATASET_LABEL\t{label}\n")
        f.write(f"COLOR\t{hex_color}\n")
        f.write(f"FIELD_LABEL\t% {label}\n")
        f.write(f"FIELD_COLORS\t{hex_color}\n")
        f.write("DATA\n")
        for name, val in sorted(data.items()):
            f.write(f"{name}\t{hex_color}\n")

# === MAIN FUNCTION ===
def build_tree_and_itol_annotations(proteome_csv_path: Path,
                                    sam_data_dir: Path,
                                    prediction_dir: Path,
                                    output_dir: Path):
    output_dir.mkdir(exist_ok=True)

    # Load proteome file
    df = pd.read_csv(proteome_csv_path)
    df["safe_label"] = df["Organism"].apply(safe_label)

    # Build tree from taxonomic lineage
    root = TreeNode()
    for _, row in df.iterrows():
        nodes = lineage_to_nodes_all(row["Taxonomic lineage"])
        root.add_path(nodes, leaf_label=row["safe_label"])
    newick = to_newick(root, is_root=True) + ";"

    # Save tree
    newick_path = output_dir / "mycobacterium_tree.nwk"
    with open(newick_path, "w") as f:
        f.write(newick)

    # Gather stats
    annotation_stats = defaultdict(lambda: {"num_sams": 0, "total_proteins": 0,
                                            "RNA": 0, "DNA": 0, "Protein": 0, "Other": 0})

    for file in sam_data_dir.glob("*.csv"):
        try:
            df_sam = pd.read_csv(file)
        except Exception:
            continue
        if "Organism" not in df_sam.columns:
            continue
        for _, row in df_sam.iterrows():
            org = safe_label(row["Organism"])
            annotation_stats[org]["num_sams"] += 1
            annotation_stats[org]["total_proteins"] += 1

    for file in prediction_dir.glob("*.csv"):
        try:
            df_pred = pd.read_csv(file)
        except Exception:
            continue
        if not {"Organism", "RNA", "DNA", "Protein", "Other"}.issubset(df_pred.columns):
            continue
        for _, row in df_pred.iterrows():
            org = safe_label(row["Organism"])
            for label in ["RNA", "DNA", "Protein", "Other"]:
                if pd.notna(row[label]) and str(row[label]) == "1":
                    annotation_stats[org][label] += 1

    # Collect percentage maps
    pct_sam_map = {}
    pct_rna_map = {}
    pct_dna_map = {}
    pct_prot_map = {}
    pct_other_map = {}

    for org, stats in annotation_stats.items():
        total = stats["total_proteins"]
        sam_total = stats["num_sams"]
        if total == 0 or sam_total == 0:
            continue
        pct_sam_map[org] = (sam_total / total) * 100
        pct_rna_map[org] = (stats["RNA"] / sam_total) * 100
        pct_dna_map[org] = (stats["DNA"] / sam_total) * 100
        pct_prot_map[org] = (stats["Protein"] / sam_total) * 100
        pct_other_map[org] = (stats["Other"] / sam_total) * 100

    # Write iTOL annotation files (Colorstrip format)
    write_itol_colorstrip(pct_sam_map, output_dir / "itol_sam_binding.txt", "SAM MTase", "#2b8cbe")
    write_itol_colorstrip(pct_rna_map, output_dir / "itol_rna_binding.txt", "RNA-binding", "#cb181d")
    write_itol_colorstrip(pct_dna_map, output_dir / "itol_dna_binding.txt", "DNA-binding", "#6a51a3")
    write_itol_colorstrip(pct_prot_map, output_dir / "itol_protein_binding.txt", "Protein-binding", "#238b45")
    write_itol_colorstrip(pct_other_map, output_dir / "itol_other_binding.txt", "Other-binding", "#d94801")

    return newick_path

build_tree_and_itol_annotations(
    proteome_csv_path=Path("../reference_files/reference_proteomes.csv"),
    sam_data_dir=Path("../reference_files"),
    prediction_dir=Path("../reference_files/prediction_files"),
    output_dir=Path("./itol_output_reference")
)