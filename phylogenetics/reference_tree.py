# =============================================================
# iTOL PIPELINE: Mycobacterium SAM Methyltransferase Tree
# =============================================================
# This script:
# 1. Builds a Newick tree from mycobacterium_proteomes.csv
# 2. Collects SAM counts and binding predictions per organism
# 3. Filters to organisms with sufficient SAM data
# 4. Prunes tree to matching organisms
# 5. Writes iTOL annotation files (DATASET_GRADIENT, DATASET_SIMPLEBAR)
# =============================================================

import pandas as pd
import re
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict
from io import StringIO

# ================= USER PATHS =================
PROTEOME_CSV = Path("../reference_files/reference_proteomes.csv")
SAM_DIR = Path("../reference_files/reference_proteomes")
PRED_DIR = Path("../reference_files/prediction_files")
OUTDIR = Path("itol_output_reference")

# Minimum SAMs required per organism to be included
MIN_SAMS = 20

# ================= HELPERS =================
def canon(name: str) -> str:
    """
    Create a canonical ID for both tree leaves and annotations.
    Ensures consistent naming across all data sources.
    """
    if not isinstance(name, str) or not name.strip():
        return "unnamed"
    # Replace whitespace with underscores
    name = re.sub(r"\s+", "_", name.strip())
    # Keep only alphanumeric, underscore, hyphen (remove brackets, parens, etc.)
    name = re.sub(r"[^A-Za-z0-9_\-]+", "", name)
    # Collapse multiple underscores
    name = re.sub(r"_+", "_", name).strip("_")
    return name if name else "unnamed"


def split_lineage(lineage) -> List[str]:
    """
    Split a taxonomic lineage string into canonical node names.
    Handles both comma and semicolon separators.
    Strips rank annotations like '(phylum)' before canonicalizing.
    """
    if not isinstance(lineage, str) or not lineage.strip():
        return []

    # Determine separator
    sep = "," if lineage.count(",") > lineage.count(";") else ";"

    parts = []
    for part in lineage.split(sep):
        part = part.strip()
        if not part:
            continue
        # Remove rank annotations like "(phylum)", "(class)", etc.
        clean = re.sub(r"\s*\([^)]*\)\s*", "", part).strip()
        if clean:
            parts.append(canon(clean))

    return parts


# ================= TREE STRUCTURE =================
class Node:
    """Simple tree node for building Newick trees."""

    def __init__(self, name: Optional[str] = None):
        self.name = name
        self.children: Dict[str, "Node"] = {}

    def add_path(self, path: List[str], leaf: str):
        """Add a path through the tree, ending at a leaf node."""
        node = self
        for p in path:
            if p not in node.children:
                node.children[p] = Node(p)
            node = node.children[p]
        # Add leaf as child of the deepest path node
        if leaf and leaf not in node.children:
            node.children[leaf] = Node(leaf)


def to_newick(node: Node, branch_len: float = 1.0) -> str:
    """
    Convert a Node tree to Newick format string.
    """
    if not node.children:
        # Leaf node
        name = node.name if node.name else ""
        return f"{name}:{branch_len}"

    # Internal node - recurse on children
    children_str = ",".join(
        to_newick(child, branch_len) for child in node.children.values()
    )

    # Root node (name=None) gets no label, internal nodes get their name
    if node.name:
        return f"({children_str}){node.name}:{branch_len}"
    else:
        return f"({children_str})"


def get_leaves(node: Node) -> List[str]:
    """Recursively collect all leaf names from a Node tree."""
    if not node.children:
        return [node.name] if node.name else []

    leaves = []
    for child in node.children.values():
        leaves.extend(get_leaves(child))
    return leaves


# ================= iTOL WRITERS =================
def write_gradient(fname: Path, label: str, values: Dict[str, float], color: str = "#084594"):
    """
    Write an iTOL DATASET_GRADIENT file.
    """
    if not values:
        print(f"  Skipping {fname.name} - no data")
        return

    vmin = min(values.values())
    vmax = max(values.values())

    with open(fname, "w") as f:
        f.write("DATASET_GRADIENT\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"DATASET_LABEL\t{label}\n")
        f.write(f"COLOR\t{color}\n")
        f.write("COLOR_MIN\t#f7fbff\n")
        f.write(f"COLOR_MAX\t{color}\n")
        f.write("DATA\n")
        for k, v in sorted(values.items()):
            f.write(f"{k}\t{v:.4f}\n")

    print(f"  Wrote {fname.name} ({len(values)} entries, range: {vmin:.3f}-{vmax:.3f})")


def write_bar(fname: Path, label: str, values: Dict[str, int], color: str = "#2171b5"):
    """
    Write an iTOL DATASET_SIMPLEBAR file.
    """
    if not values:
        print(f"  Skipping {fname.name} - no data")
        return

    with open(fname, "w") as f:
        f.write("DATASET_SIMPLEBAR\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"DATASET_LABEL\t{label}\n")
        f.write(f"COLOR\t{color}\n")
        f.write("DATA\n")
        for k, v in sorted(values.items()):
            f.write(f"{k}\t{v}\n")

    print(f"  Wrote {fname.name} ({len(values)} entries)")


# ================= MAIN PIPELINE =================
def main():
    # --- Validate inputs ---
    if not PROTEOME_CSV.exists():
        raise FileNotFoundError(f"Proteome CSV not found: {PROTEOME_CSV}")
    if not SAM_DIR.exists():
        raise FileNotFoundError(f"SAM data directory not found: {SAM_DIR}")
    if not PRED_DIR.exists():
        raise FileNotFoundError(f"Prediction directory not found: {PRED_DIR}")

    OUTDIR.mkdir(exist_ok=True)

    # --- Load proteome metadata ---
    print("Loading proteome metadata...")
    df = pd.read_csv(PROTEOME_CSV)

    required_cols = {"Organism", "Taxonomic lineage"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Proteome CSV missing required columns: {required_cols - set(df.columns)}")

    df["tree_id"] = df["Organism"].apply(canon)
    print(f"  Loaded {len(df)} organisms from {PROTEOME_CSV.name}")

    # --- Build tree from taxonomic lineages ---
    print("\nBuilding taxonomic tree...")
    root = Node(None)  # Root has no name

    for _, row in df.iterrows():
        lineage_nodes = split_lineage(row["Taxonomic lineage"])
        tree_id = row["tree_id"]

        if lineage_nodes:
            root.add_path(lineage_nodes, tree_id)
        else:
            # No lineage - add directly under root
            root.add_path([], tree_id)

    original_leaves = get_leaves(root)
    print(f"  Tree has {len(original_leaves)} leaves")

    # --- Collect SAM and prediction statistics ---
    print("\nCollecting SAM statistics...")
    stats = defaultdict(lambda: {
        "total": 0,
        "sams": 0,
        "RNA": 0,
        "DNA": 0,
        "Protein": 0,
        "Other": 0
    })

    # Process SAM data files
    sam_files = list(SAM_DIR.glob("*.csv"))
    print(f"  Found {len(sam_files)} CSV files in {SAM_DIR}")

    for f in sam_files:
        try:
            d = pd.read_csv(f, low_memory=False)
        except Exception as e:
            print(f"    Could not read {f.name}: {e}")
            continue

        if "Organism" not in d.columns:
            continue

        # Get organism name from first row
        org_raw = d["Organism"].iloc[0] if len(d) > 0 else None
        if pd.isna(org_raw):
            continue

        org = canon(str(org_raw))
        stats[org]["total"] += len(d)

        # Count SAMs (proteins with IPR029063 InterPro domain)
        if "InterPro" in d.columns:
            sam_count = d["InterPro"].fillna("").str.contains("IPR029063", regex=False).sum()
            stats[org]["sams"] += sam_count

    print(f"  Collected stats for {len(stats)} organisms")

    # Process prediction files
    print("\nCollecting binding predictions...")
    pred_files = list(PRED_DIR.glob("*.csv"))
    print(f"  Found {len(pred_files)} CSV files in {PRED_DIR}")

    for f in pred_files:
        try:
            d = pd.read_csv(f, low_memory=False)
        except Exception as e:
            print(f"    Could not read {f.name}: {e}")
            continue

        required = {"Organism", "RNA", "DNA", "Protein"}
        if not required.issubset(d.columns):
            continue

        org_raw = d["Organism"].iloc[0] if len(d) > 0 else None
        if pd.isna(org_raw):
            continue

        org = canon(str(org_raw))

        # Accumulate binding predictions (threshold > 0.5)
        for col in ["RNA", "DNA", "Protein", "Other"]:
            if col in d.columns:
                # Convert to numeric, count values > 0.5
                values = pd.to_numeric(d[col], errors="coerce")
                stats[org][col] += (values > 0.5).sum()

    # --- Filter organisms by minimum SAM count ---
    print(f"\nFiltering organisms (minimum {MIN_SAMS} SAMs)...")

    valid_organisms = set()
    pct_sams = {}
    n_sams = {}
    pct_RNA = {}
    pct_DNA = {}
    pct_PROT = {}
    pct_OTHER = {}

    for org, s in stats.items():
        if s["sams"] >= MIN_SAMS and s["total"] > 0:
            valid_organisms.add(org)

            pct_sams[org] = s["sams"] / s["total"]
            n_sams[org] = s["sams"]

            # Calculate binding percentages (relative to SAM count)
            sam_count = s["sams"]
            pct_RNA[org] = s["RNA"] / sam_count if sam_count > 0 else 0.0
            pct_DNA[org] = s["DNA"] / sam_count if sam_count > 0 else 0.0
            pct_PROT[org] = s["Protein"] / sam_count if sam_count > 0 else 0.0
            pct_OTHER[org] = s["Other"] / sam_count if sam_count > 0 else 0.0

    print(f"  {len(valid_organisms)} organisms passed filter (of {len(stats)} total)")

    if not valid_organisms:
        print("\nNo organisms passed the filter. Check your data or lower MIN_SAMS.")
        return

    # --- Prune tree to valid organisms ---
    print("\nPruning tree to valid organisms...")

    def prune_tree(node: Node, valid: set) -> Optional[Node]:
        """Recursively prune tree, keeping only branches leading to valid leaves."""
        if not node.children:
            # Leaf node - keep only if valid
            return node if node.name in valid else None

        # Internal node - prune children
        new_children = {}
        for name, child in node.children.items():
            pruned = prune_tree(child, valid)
            if pruned is not None:
                new_children[name] = pruned

        if not new_children:
            return None

        # Create pruned node
        pruned_node = Node(node.name)
        pruned_node.children = new_children
        return pruned_node

    pruned_root = prune_tree(root, valid_organisms)

    if pruned_root is None:
        print("  Tree is empty after pruning!")
        return

    pruned_leaves = get_leaves(pruned_root)
    print(f"  Pruned tree has {len(pruned_leaves)} leaves")

    # --- Write Newick tree ---
    nwk_path = OUTDIR / "reference_tree.nwk"
    newick_str = to_newick(pruned_root) + ";"
    nwk_path.write_text(newick_str)
    print(f"  Wrote tree to {nwk_path}")

    # --- Filter annotation maps to match tree leaves ---
    tree_leaf_set = set(pruned_leaves)

    pct_sams = {k: v for k, v in pct_sams.items() if k in tree_leaf_set}
    n_sams = {k: v for k, v in n_sams.items() if k in tree_leaf_set}
    pct_RNA = {k: v for k, v in pct_RNA.items() if k in tree_leaf_set}
    pct_DNA = {k: v for k, v in pct_DNA.items() if k in tree_leaf_set}
    pct_PROT = {k: v for k, v in pct_PROT.items() if k in tree_leaf_set}
    pct_OTHER = {k: v for k, v in pct_OTHER.items() if k in tree_leaf_set}

    # --- Write iTOL annotation files ---
    print("\nWriting iTOL annotation files...")

    write_gradient(OUTDIR / "itol_pct_sams.txt", "% SAMs in proteome", pct_sams, "#2171b5")
    write_bar(OUTDIR / "itol_n_sams.txt", "# SAM MTases", n_sams, "#2b8cbe")
    write_gradient(OUTDIR / "itol_pct_rna.txt", "% RNA-binding SAMs", pct_RNA, "#cb181d")
    write_gradient(OUTDIR / "itol_pct_dna.txt", "% DNA-binding SAMs", pct_DNA, "#6a51a3")
    write_gradient(OUTDIR / "itol_pct_protein.txt", "% Protein-binding SAMs", pct_PROT, "#238b45")
    write_gradient(OUTDIR / "itol_pct_other.txt", "% Other SAMs", pct_OTHER, "#e6550d")

    # --- Final report ---
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"Input organisms:        {len(df)}")
    print(f"Organisms with data:    {len(stats)}")
    print(f"Passed SAM filter:      {len(valid_organisms)}")
    print(f"Final tree leaves:      {len(pruned_leaves)}")
    print(f"Output directory:       {OUTDIR.resolve()}")
    print("=" * 50)

    # Show per-organism stats
    print("\nPer-organism statistics:")
    print(f"{'Organism':<40} {'SAMs':>6} {'%RNA':>7} {'%DNA':>7} {'%Prot':>7} {'%Other':>7}")
    print("-" * 80)
    for org in sorted(pruned_leaves):
        if org in n_sams:
            print(f"{org:<40} {n_sams[org]:>6} "
                  f"{pct_RNA.get(org, 0)*100:>6.1f}% "
                  f"{pct_DNA.get(org, 0)*100:>6.1f}% "
                  f"{pct_PROT.get(org, 0)*100:>6.1f}% "
                  f"{pct_OTHER.get(org, 0)*100:>6.1f}%")


if __name__ == "__main__":
    main()
