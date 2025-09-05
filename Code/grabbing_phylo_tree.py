# import os, re
# import pandas as pd
# from ete3 import NCBITaxa

# # ============ CONFIG ============
# FOLDER = "prediction_files"
# OUT_TREE = "species_tree.nwk"
# OUT_UNRESOLVED = "ncbi_not_found.csv"
# OUT_AUDIT = "ncbi_resolution_audit.csv"
# BATCH = 200  # batch size for name_translator

# # ============ REGEX / CLEANING ============
# # Expand culture collection tokens and qualifiers
# COLLECTION_TOKENS = (
#     r"(?:ATCC|DSM|JCM|NBRC|NRRL|NCTC|CCUG|CIP|BCRC|KCTC|CGMCC|MCCC|PCC|CBS|FGSC|IAM|"
#     r"CCALA|CCMP|CCAP|MBIC|RCEF|IBT|ARh|SAG|RIKEN|NCIMB|UWO|NRZ|IFO|NBIMCC|NIES|KMM|"
#     r"KACC|LMG)"
#     r"|(?:strain|isolate|serotype|serogroup|biovar|biotype|lineage|sublineage|clade)"
# )
# QUALIFIER_TOKENS = r"(?:subsp\.|ssp\.|var\.|forma|subspecies|pv\.)"

# RE_PARENS_ALL = re.compile(r"\s*\([^)]*\)")
# RE_QUOTES = re.compile(r"[‘’“”'\"`]+")
# RE_AFTER_EQUAL = re.compile(r"\s*=\s*.*$")                   # cut at " = ..."
# RE_COLLECTION = re.compile(rf"\b{COLLECTION_TOKENS}\b.*$", re.IGNORECASE)
# RE_SP_CODE = re.compile(r"\bsp\.?\b.*$", re.IGNORECASE)     # drop "sp. XXX"
# RE_FORM_PV = re.compile(rf"\b(?:f\.?\s*sp\.?|f\.|pv\.)\b.*$", re.IGNORECASE)  # drop forma/pathovar tail
# RE_SEROVAR = re.compile(r"\bserovar\b.*$", re.IGNORECASE)
# RE_SEROGROUP = re.compile(r"\bserogroup\b.*$", re.IGNORECASE)
# RE_HYBRID = re.compile(r"\b[x×]\b")                         # hybrids "A x B"
# RE_MULTI_WS = re.compile(r"\s+")
# RE_DANGLING = re.compile(r"^[\W_]+|[\W_]+$")
# RE_TRAIL_PARENS = re.compile(r"[()]+$")
# RE_SQUOTE_DASH = re.compile(r"[’´`]+")

# # Trailing tokens that look like lab/isolate codes even without a keyword
# # e.g., "FD-317", "FP15055", "FSL", "1063", "ZS01", "RBIITD", "JGI-2019a"
# RE_TRAILING_CODE = re.compile(
#     r"(?:\s(?:[A-Z]{1,6}(?:[-_]*\d+[A-Z]*)|\d{3,}|[A-Z]{2,}[-_]\d{2,}))$"
# )

# # Optional: strip multiple trailing codes iteratively
# def strip_trailing_codes(name: str) -> str:
#     prev = None
#     s = name
#     # Strip up to 3 times to be safe (e.g., "FD-317 2021")
#     for _ in range(3):
#         if prev == s:
#             break
#         prev = s
#         s = re.sub(r"\s+" + RE_TRAILING_CODE.pattern, "", s).strip()
#     return s

# def normalize_spaces(s: str) -> str:
#     return s.replace("\xa0", " ").replace("\u2009", " ").replace("\u2002", " ").replace("\u2003", " ")

# def tidy_name(raw: str) -> str:
#     """
#     Aggressive cleaning to produce a scientific name candidate:
#     - remove quotes, parentheticals, equivalences, culture codes, form/pathovar tails, serovar/group
#     - collapse sp. isolates to genus
#     - handle hybrids "A x B" by keeping genus
#     - bracketed genera -> unbracketed
#     """
#     if not isinstance(raw, str) or not raw.strip():
#         return ""
#     name = normalize_spaces(raw.strip())
#     name = name.replace("_", " ")
#     name = RE_QUOTES.sub("", name)
#     name = RE_SQUOTE_DASH.sub("", name)

#     # keep only left side of equivalences
#     name = RE_AFTER_EQUAL.sub("", name)
#     # remove all parentheticals (common names, synonyms, strain notes)
#     name = RE_PARENS_ALL.sub("", name)

#     # standardize bracketed genus: [Genus] species -> Genus species
#     name = name.replace("[", "").replace("]", "")

#     # drop form/pathovar tails & serovar/serogroup tails
#     name = RE_FORM_PV.sub("", name)
#     name = RE_SEROVAR.sub("", name)
#     name = RE_SEROGROUP.sub("", name)

#     # cut culture/strain tokens onwards
#     name = RE_COLLECTION.sub("", name)

#     # cut after "sp." and any trailing isolate code
#     name = RE_SP_CODE.sub("", name)

#     # hybrids: if "A x B", reduce to genus (first token)
#     if RE_HYBRID.search(name):
#         toks = name.split()
#         if toks:
#             name = toks[0]  # Genus level fallback

#     # compress whitespace, trim punctuation/parens
#     name = RE_MULTI_WS.sub(" ", name).strip()
#     name = RE_TRAIL_PARENS.sub("", name)
#     name = RE_DANGLING.sub("", name).strip()

#     # Handle qualifiers (subsp./var./forma): keep up to 4 tokens (Genus species qualifier epithet), else 2–3 tokens
#     toks = name.split()
#     if len(toks) >= 3 and toks[2].lower() in {"subsp.", "ssp.", "var.", "forma", "subspecies", "pv."}:
#         name = " ".join(toks[:4]) if len(toks) >= 4 else " ".join(toks[:3])
#     else:
#         if len(toks) > 3:
#             name = " ".join(toks[:3])

#     # "spp." or "sp." endings -> keep genus only
#     low = name.lower()
#     if low.endswith((" sp", " sp.", " spp", " spp.")):
#         name = name.split()[0]

#     # remove trailing isolate-like codes with no keyword
#     name = strip_trailing_codes(name)
#     # recollapse whitespace
#     name = RE_MULTI_WS.sub(" ", name).strip()

#     return name

# # Optional targeted synonym/normalization hooks (kept conservative)
# def synonym_normalize(name: str) -> str:
#     low = name.lower()

#     # General phytoplasma: if the word appears anywhere, map to genus-level
#     if "phytoplasma" in low:
#         return "Candidatus Phytoplasma"

#     # Common truncated host phrases ending with 'aster' from 'aster yellows phytoplasma'
#     if low.endswith(" aster") or "aster yellows" in low:
#         return "Candidatus Phytoplasma asteris"

#     # Endosymbiont phrasing: "Cardinium endosymbiont of ..." -> "Candidatus Cardinium"
#     if "endosymbiont" in low:
#         first = name.split()[0]
#         if not first.lower().startswith("candidatus"):
#             return f"Candidatus {first}"
#         return first

#     # Ensure Candidatus is properly capitalized if present
#     if low.startswith("candidatus "):
#         return "Candidatus " + name.split(" ", 1)[1]

#     return name

# def best_taxid_from_search(ncbi: NCBITaxa, q: str):
#     try:
#         hits = ncbi.search_taxa(q)
#     except Exception:
#         hits = []
#     if not hits:
#         return None
#     # exact case-insensitive first
#     for h in hits:
#         if h.get("name", "").lower() == q.lower():
#             return h["taxid"]
#     # prefer species/subspecies
#     for h in hits:
#         if h.get("rank") in {"species", "subspecies"}:
#             return h["taxid"]
#     # then genus
#     for h in hits:
#         if h.get("rank") == "genus":
#             return h["taxid"]
#     return hits[0]["taxid"]

# def main():
#     ncbi = NCBITaxa()
#     # Uncomment to refresh local taxonomy (one-time maintenance):
#     # ncbi.update_taxonomy_database()

#     # 1) Gather raw names
#     raw_names = set()
#     for fname in os.listdir(FOLDER):
#         if not fname.endswith(".csv"):
#             continue
#         fpath = os.path.join(FOLDER, fname)
#         try:
#             df = pd.read_csv(fpath, dtype=str, encoding_errors="ignore", on_bad_lines="skip")
#         except Exception as e:
#             print(f"⚠️ Could not read {fname}: {e}")
#             continue
#         col = next((c for c in ["Organism", "organism", "Species", "species"] if c in df.columns), None)
#         if not col:
#             print(f"⚠️ No 'Organism' column in {fname}")
#             continue
#         raw_names.update(df[col].dropna().astype(str).tolist())

#     print(f"Collected {len(raw_names)} raw names.")

#     # 2) Clean + normalize
#     pairs = []
#     for r in raw_names:
#         c = tidy_name(r)
#         c = synonym_normalize(c) if c else c
#         if c:
#             pairs.append((r, c))
#         else:
#             pairs.append((r, ""))

#     # 3) Batch exact translate
#     cleaned_unique = sorted({c for _, c in pairs if c})
#     name2tax = {}
#     for i in range(0, len(cleaned_unique), BATCH):
#         batch = cleaned_unique[i:i+BATCH]
#         try:
#             tr = ncbi.get_name_translator(batch)  # name -> [taxid,...]
#         except Exception:
#             tr = {}
#         for n in batch:
#             if n in tr and tr[n]:
#                 name2tax[n] = tr[n][0]

#     # 4) Resolve with fallbacks + collect audit
#     rows = []
#     unresolved = []
#     taxids = set()

#     for raw, cn in pairs:
#         matched_name = None
#         rank = None
#         tid = None

#         if cn:
#             tid = name2tax.get(cn)

#             if tid is None:
#                 # try species/subspecies; then first two tokens; then genus
#                 tid = best_taxid_from_search(ncbi, cn)
#                 if tid is None:
#                     toks = cn.split()
#                     if len(toks) >= 2:
#                         tid = best_taxid_from_search(ncbi, " ".join(toks[:2]))
#                 if tid is None and cn:
#                     tid = best_taxid_from_search(ncbi, cn.split()[0])

#         # if tid is None and cn:
#         #     # final attempt: if qualifiers existed earlier, try base species (first two tokens)
#         #     toks = cn.split()
#         #     if len(toks) >= 2:
#         #         tid = best_taxid_from_search(ncbi, " ".join(toks[:2]))

#         if tid is None:
#             unresolved.append((raw, cn))
#             rows.append([raw, cn, None, None, None])
#         else:
#             # fetch scientific name & rank for audit
#             try:
#                 name = ncbi.get_taxid_translator([tid]).get(tid, None)
#             except Exception:
#                 name = None
#             try:
#                 ranks = ncbi.get_rank([tid])
#                 rank = ranks.get(tid, None)
#             except Exception:
#                 rank = None
#             matched_name = name
#             taxids.add(tid)
#             rows.append([raw, cn, tid, matched_name, rank])

#     # 5) Write audits
#     pd.DataFrame(rows, columns=["raw_name", "cleaned_query", "taxid", "matched_name", "rank"]).to_csv(OUT_AUDIT, index=False)
#     if unresolved:
#         pd.DataFrame(unresolved, columns=["raw_name", "cleaned_query"]).to_csv(OUT_UNRESOLVED, index=False)
#         print(f"Unresolved: {len(unresolved)} → {OUT_UNRESOLVED}")
#     print(f"Resolved: {len(taxids)} unique taxids. Audit → {OUT_AUDIT}")

#     # 6) Build tree
#     if taxids:
#         tree = ncbi.get_topology(sorted(taxids))
#         tree.write(format=1, outfile=OUT_TREE)
#         print(f"🌳 Tree written → {OUT_TREE} (n={len(taxids)} taxa)")
#     else:
#         print("❌ No valid taxids to build a tree.")

# if __name__ == "__main__":
#     main()

# #!/usr/bin/env python3
# import os
# import re
# import csv
# import gzip
# import io
# import argparse
# from pathlib import Path
# from collections import Counter

# import pandas as pd

# # ==================== CONFIG (defaults; can override via CLI) ====================
# FOLDER = "prediction_files"   # where your .csvs live
# EXCEL_GZ = "Files/proteomes_AND_proteome_type_1_2025_08_28.xlsx.gz"  # gzip'd Excel with Organism & Taxonomic lineage

# OUT_NEWICK = "taxonomy_tree.nwk"
# OUT_USED = "used_lineages.csv"              # organism,lineage actually used
# OUT_UNMATCHED = "unmatched_organisms.csv"   # organisms we couldn't match to a lineage

# DROP_SUPERROOTS = {"cellular organisms"}    # drop this super-root if present
# USE_SPECIES_AS_LEAVES = True                # attach organism as leaf under last taxon node
# BRANCH_LENGTH = 1                           # constant branch length for Newick
# EXCEL_SHEET = None                          # None = first sheet

# # ================== HELPERS =====================
# def vprint(verbose, *args):
#     if verbose:
#         print(*args)

# def norm_org_key(s: str) -> str:
#     """Normalize organism string for matching."""
#     if not isinstance(s, str):
#         return ""
#     s = s.replace("_", " ").replace("\xa0", " ").replace("’", "'").strip()
#     s = re.sub(r"\s+", " ", s)
#     s = s.strip('"\''"“”‘’")
#     return s.lower()

# def split_lineage(s: str, sep_guess=None):
#     if not isinstance(s, str):
#         return []
#     s = s.strip()
#     if not s:
#         return []
#     sep = sep_guess
#     if sep is None:
#         c_comma = s.count(",")
#         c_semi  = s.count(";")
#         sep = "," if c_comma >= c_semi else ";"
#     parts = [p.strip() for p in s.split(sep)]
#     return [p for p in parts if p]

# def strip_rank(tag: str):
#     if not isinstance(tag, str):
#         return ""
#     m = re.match(r"^(.*?)(?:\s*\([^)]*\))\s*$", tag)
#     name = m.group(1) if m else tag
#     name = re.sub(r"\s+", " ", name).strip()
#     return name

# def clean_taxon_name(name: str):
#     name = strip_rank(name)
#     name = name.replace('"', "").replace("'", "").replace("\xa0", " ").strip()
#     return name

# def lineage_to_nodes(lineage_str: str):
#     parts = split_lineage(lineage_str)
#     parts = [clean_taxon_name(p) for p in parts]
#     parts = [p for p in parts if p]
#     parts = [p for p in parts if p not in DROP_SUPERROOTS]
#     dedup = []
#     for p in parts:
#         if not dedup or dedup[-1] != p:
#             dedup.append(p)
#     return dedup

# class TreeNode:
#     __slots__ = ("name", "children", "is_leaf_tag")
#     def __init__(self, name=None):
#         self.name = name
#         self.children = {}  # name -> TreeNode
#         self.is_leaf_tag = False

#     def add_path(self, path, leaf_label=None):
#         node = self
#         for n in path:
#             if n not in node.children:
#                 node.children[n] = TreeNode(n)
#             node = node.children[n]
#         if leaf_label is not None:
#             if leaf_label not in node.children:
#                 node.children[leaf_label] = TreeNode(leaf_label)
#             node.children[leaf_label].is_leaf_tag = True

# def escape_newick_label(label: str):
#     if label is None:
#         return ""
#     label = label.replace(" ", "_")
#     label = re.sub(r"[,:;()]", "_", label)
#     return label

# def to_newick(node: TreeNode):
#     if not node.children:
#         label = escape_newick_label(node.name) if node.name is not None else ""
#         return f"{label}:{BRANCH_LENGTH}"
#     subs = ",".join(to_newick(child) for child in node.children.values())
#     label = "" if node.name is None else escape_newick_label(node.name)
#     return f"({subs}){label}:{BRANCH_LENGTH}"

# # ================== PIPELINE FUNCTIONS ====================
# def collect_unique_organisms(folder: str, verbose=False) -> set:
#     orgs = set()
#     folder_path = Path(folder)
#     if not folder_path.exists():
#         raise FileNotFoundError(f"Folder not found: {folder}")
#     files = sorted([f for f in folder_path.iterdir() if f.suffix.lower() == ".csv"])
#     if not files:
#         print(f"⚠️ No CSVs found under {folder}")
#     vprint(verbose, f"Found {len(files)} CSV files in {folder}")
#     for fpath in files:
#         try:
#             df = pd.read_csv(fpath, dtype=str, encoding_errors="ignore", on_bad_lines="skip")
#         except Exception as e:
#             print(f"⚠️ Could not read {fpath.name}: {e}")
#             continue
#         col = next((c for c in ["Organism", "organism", "Species", "species"] if c in df.columns), None)
#         if not col:
#             print(f"⚠️ No 'Organism' column in {fpath.name} (columns: {list(df.columns)[:8]}...)")
#             continue
#         vals = df[col].dropna().astype(str).tolist()
#         before = len(orgs)
#         for v in vals:
#             v2 = v.strip()
#             if v2:
#                 orgs.add(v2)
#         added = len(orgs) - before
#         vprint(verbose, f"  {fpath.name}: +{added} new organisms (total {len(orgs)})")
#     return orgs

# def load_lineage_table_from_excel_gz(path: str, excel_sheet=None, verbose=False) -> pd.DataFrame:
#     p = Path(path)
#     if not p.exists():
#         raise FileNotFoundError(f"Lineage file not found: {path}")

#     with gzip.open(p, "rb") as fh:
#         data = fh.read()
#     bio = io.BytesIO(data)

#     # read first sheet or specified sheet
#     # If excel_sheet is None, use the FIRST sheet (index 0), not "all sheets"
#     sheet_kw = 0 if excel_sheet is None else excel_sheet

#     df = pd.read_excel(bio, sheet_name=sheet_kw, dtype=str, engine="openpyxl")

#     # If a dict slipped through (e.g., user passed --sheet None), pick the first sheet
#     if isinstance(df, dict):
#         # preserve order; pick the first by insertion order
#         first_name = next(iter(df.keys()))
#         if verbose:
#             print(f"ℹ️ Multiple sheets detected; using first: {first_name}")
#         df = df[first_name]

#     cols_lower = [c.lower() for c in df.columns]
#     if verbose:
#         print(f"Lineage Excel columns: {cols_lower}")

#     # normalize columns (case-insensitive)
#     cols = {c.lower(): c for c in df.columns}
#     org_col = cols.get("organism")
#     lin_col = cols.get("taxonomic lineage") or cols.get("taxonomic_lineage") or cols.get("lineage")
#     if not org_col or not lin_col:
#         raise ValueError("Excel must have columns: 'Organism' and 'Taxonomic lineage' (case-insensitive)")

#     sub = df[[org_col, lin_col]].copy()
#     sub.columns = ["Organism", "Taxonomic lineage"]
#     sub["__key__"] = sub["Organism"].astype(str).map(norm_org_key)
#     # Drop duplicate organisms in the lineage table (keep first)
#     before = len(sub)
#     sub = sub.drop_duplicates(subset="__key__", keep="first")
#     vprint(verbose, f"Lineage rows: {before} -> {len(sub)} after de-duplication by organism key")
#     return sub

# def match_orgs_to_lineages(unique_orgs: set, lineage_df: pd.DataFrame, verbose=False):
#     used = []
#     unmatched = []
#     key2lineage = dict(zip(lineage_df["__key__"], lineage_df["Taxonomic lineage"]))

#     for org in sorted(unique_orgs):
#         key = norm_org_key(org)
#         lin = key2lineage.get(key)
#         if lin is None:
#             # light fallback: underscores->spaces already handled in norm; try as-is again
#             lin = key2lineage.get(norm_org_key(org.replace("_", " ")))
#         if lin is None:
#             unmatched.append(org)
#         else:
#             used.append((org, lin))
#     vprint(verbose, f"Matched {len(used)} / {len(unique_orgs)} organisms to lineages")
#     return used, unmatched

# def build_tree_from_used(used_pairs, verbose=False):
#     root = TreeNode(name=None)
#     dup_counter = Counter()
#     attached = 0
#     for org, lineage in used_pairs:
#         path = lineage_to_nodes(lineage)
#         if not path:
#             continue
#         if USE_SPECIES_AS_LEAVES:
#             leaf = org
#             dup_counter[leaf] += 1
#             if dup_counter[leaf] > 1:
#                 leaf = f"{leaf}__{dup_counter[org]}"
#             root.add_path(path, leaf_label=leaf)
#         else:
#             root.add_path(path)
#         attached += 1
#     vprint(verbose, f"Attached {attached} leaves to tree")
#     newick = to_newick(root)
#     if not newick.endswith(";"):
#         newick += ";"
#     return newick

# # ================== MAIN ====================
# def main():
#     ap = argparse.ArgumentParser(description="Build Newick tree from lineage table + prediction_files organisms.")
#     ap.add_argument("--folder", default=FOLDER, help="Path to folder with prediction CSVs")
#     ap.add_argument("--excel", default=EXCEL_GZ, help="Path to gzipped Excel with Organism + Taxonomic lineage")
#     ap.add_argument("--sheet", default=EXCEL_SHEET, help="Excel sheet name (default: first)")
#     ap.add_argument("--out-newick", default=OUT_NEWICK)
#     ap.add_argument("--out-used", default=OUT_USED)
#     ap.add_argument("--out-unmatched", default=OUT_UNMATCHED)
#     ap.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
#     args = ap.parse_args()

#     print("=== Lineage → Newick Builder ===")
#     print(f"CSV folder:       {args.folder}")
#     print(f"Lineage Excel.gz: {args.excel}")
#     print(f"Sheet:            {args.sheet or '(first)'}")

#     # 1) gather organisms
#     unique_orgs = collect_unique_organisms(args.folder, verbose=args.verbose)
#     print(f"✅ Collected {len(unique_orgs)} unique organisms")

#     # 2) load lineage table
#     lineage_df = load_lineage_table_from_excel_gz(args.excel, excel_sheet=args.sheet, verbose=args.verbose)
#     print(f"📚 Loaded lineage table with {len(lineage_df)} unique organisms (by normalized key)")

#     # 3) match
#     used, unmatched = match_orgs_to_lineages(unique_orgs, lineage_df, verbose=args.verbose)
#     print(f"🔗 Matched {len(used)} | ❌ Unmatched {len(unmatched)}")

#     # 4) write audits
#     if used:
#         pd.DataFrame(used, columns=["organism", "lineage"]).to_csv(args.out_used, index=False)
#         print(f"📝 Wrote used lineages → {args.out_used}")
#     if unmatched:
#         pd.DataFrame({"organism": unmatched}).to_csv(args.out_unmatched, index=False)
#         print(f"🧭 Wrote unmatched organisms → {args.out_unmatched}")

#     # 5) build tree
#     if used:
#         newick = build_tree_from_used(used, verbose=args.verbose)
#         Path(args.out_newick).write_text(newick, encoding="utf-8")
#         print(f"🌳 Wrote Newick → {args.out_newick}")
#     else:
#         print("❌ No matched organisms; Newick not created.")

# if __name__ == "__main__":
    # main()

#!/usr/bin/env python3
import os
import re
import argparse
from pathlib import Path
from collections import Counter
import pandas as pd

# ==================== CONFIG (can override via CLI) ====================
DROP_SUPERROOTS = {"cellular organisms"}   # remove this super-root if present
USE_SPECIES_AS_LEAVES = True               # attach organism as leaf under last taxon node
BRANCH_LENGTH = 1                          # constant branch length for Newick
OUT_NEWICK_DEFAULT = "taxonomy_tree.nwk"

# ==================== HELPERS (tree + lineage parsing) ====================
def split_lineage(s: str, sep_guess=None):
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

def strip_rank(tag: str):
    """Remove trailing ' (rank)' if present, keep the taxon name only."""
    if not isinstance(tag, str):
        return ""
    m = re.match(r"^(.*?)(?:\s*\([^)]*\))\s*$", tag)
    name = m.group(1) if m else tag
    name = re.sub(r"\s+", " ", name).strip()
    return name

def clean_taxon_name(name: str):
    """Basic cleaning of a taxon node label."""
    name = strip_rank(name)
    name = name.replace('"', "").replace("'", "").replace("\xa0", " ").strip()
    return name

def lineage_to_nodes(lineage_str: str):
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

def escape_newick_label(label: str):
    """Escape characters that Newick parsers dislike."""
    if label is None:
        return ""
    label = label.replace(" ", "_")
    label = re.sub(r"[,:;()]", "_", label)
    return label

def to_newick(node: TreeNode):
    """Recursively render Newick with constant branch lengths."""
    if not node.children:
        label = escape_newick_label(node.name) if node.name is not None else ""
        return f"{label}:{BRANCH_LENGTH}"
    subs = ",".join(to_newick(child) for child in node.children.values())
    label = "" if node.name is None else escape_newick_label(node.name)
    return f"({subs}){label}:{BRANCH_LENGTH}"

def build_tree_from_used(used_pairs, verbose=False):
    """
    used_pairs: list of (organism, lineage_string)
    """
    root = TreeNode(name=None)
    dup_counter = Counter()
    attached = 0
    for org, lineage in used_pairs:
        path = lineage_to_nodes(lineage)
        if not path:
            continue
        if USE_SPECIES_AS_LEAVES:
            leaf = org
            dup_counter[leaf] += 1
            if dup_counter[leaf] > 1:
                leaf = f"{leaf}__{dup_counter[org]}"
            root.add_path(path, leaf_label=leaf)
        else:
            root.add_path(path)
        attached += 1
    if verbose:
        print(f"🌿 Attached {attached} leaves to tree")
    newick = to_newick(root)
    if not newick.endswith(";"):
        newick += ";"
    return newick

def escape_newick_label_for_check(label: str) -> str:
    s = str(label).replace(" ", "_")
    s = re.sub(r"[,:;()]", "_", s)
    return s

# ==================== MAIN PIPELINE ====================
def main():
    ap = argparse.ArgumentParser(
        description="Map species to lineages via UniProt (score-1 then score-4) and build a Newick tree."
    )
    ap.add_argument("--pred-dir", default="prediction_files",
                    help="Folder with per-protein prediction CSVs (must have columns: Organism, Entry)")
    ap.add_argument("--score1", default="Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03.tsv.gz",
                    help="Gzipped TSV with columns: Entry, Taxonomic lineage (preferred)")
    ap.add_argument("--score4", default="Files/uniprotkb_IPR029063_AND_annotation_scor_2025_09_03 (1).tsv.gz",
                    help="Gzipped TSV with columns: Entry, Taxonomic lineage (fallback)")
    ap.add_argument("--tree-check", default=None,
                    help="Optional .nwk to check which species are present (leaf labels) — purely for reporting")
    ap.add_argument("--out-used", default="used_lineages_from_score14.csv",
                    help="Output CSV: species, lineage, source, n_entries_matched_1, n_entries_matched_4")
    ap.add_argument("--out-no-score14", default="species_without_score14_protein.csv",
                    help="Output CSV: species with no Entry in score-1 or score-4")
    ap.add_argument("--out-not-in-tree", default="species_not_in_tree.csv",
                    help="Output CSV: matched species not present as leaves in the provided --tree-check")
    ap.add_argument("--out-newick", default=OUT_NEWICK_DEFAULT,
                    help="Output Newick tree filename")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    pred_dir = Path(args.pred_dir)
    if not pred_dir.exists():
        raise FileNotFoundError(f"Prediction folder not found: {pred_dir}")

    # --- Load score-1 table ---
    print(f"📥 Loading score-1 table: {args.score1}")
    score1 = pd.read_csv(args.score1, sep="\t", compression="gzip", dtype=str, low_memory=False)
    needed1 = {"Entry", "Taxonomic lineage"}
    if not needed1.issubset(score1.columns):
        raise ValueError(f"{args.score1} must contain columns: {needed1} (found {list(score1.columns)[:10]} ...)")
    score1 = score1[["Entry", "Taxonomic lineage"]].dropna()
    score1["Entry"] = score1["Entry"].astype(str).str.strip()
    entry2lin1 = dict(zip(score1["Entry"], score1["Taxonomic lineage"]))
    print(f"✅ Loaded {len(entry2lin1):,} score-1 entries with lineage")

    # --- Load score-4 table (fallback) ---
    print(f"📥 Loading score-4 table: {args.score4}")
    score4 = pd.read_csv(args.score4, sep="\t", compression="gzip", dtype=str, low_memory=False)
    needed4 = {"Entry", "Taxonomic lineage"}
    if not needed4.issubset(score4.columns):
        raise ValueError(f"{args.score4} must contain columns: {needed4} (found {list(score4.columns)[:10]} ...)")
    score4 = score4[["Entry", "Taxonomic lineage"]].dropna()
    score4["Entry"] = score4["Entry"].astype(str).str.strip()
    entry2lin4 = dict(zip(score4["Entry"], score4["Taxonomic lineage"]))
    print(f"✅ Loaded {len(entry2lin4):,} score-4 entries with lineage")

    # --- Scan prediction files and collect per-species Entry IDs ---
    species_to_entries = {}
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
            if args.verbose:
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
    used_rows = []               # (species, lineage, source, n1, n4)
    no_score14_species = []      # species with zero matching Entries in score-1 or score-4
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
    not_in_tree = []
    if args.tree_check:
        tree_path = Path(args.tree_check)
        if not tree_path.exists():
            print(f"  ⚠️ Tree not found at {tree_path} — skipping tree comparison")
        else:
            nwk = tree_path.read_text(encoding="utf-8", errors="ignore")
            leaf_labels = set(re.findall(r"([^\s,():]+):", nwk))
            for sp, _, _, _, _ in used_rows:
                if escape_newick_label_for_check(sp) not in leaf_labels:
                    not_in_tree.append(sp)
            print(f"🌳 Matched species NOT present as leaves in provided tree: {len(not_in_tree):,}")

    # --- Write audits ---
    if used_rows:
        pd.DataFrame(
            used_rows,
            columns=["species", "lineage", "source", "n_entries_matched_1", "n_entries_matched_4"]
        ).to_csv(args.out_used, index=False)
        print(f"📝 Wrote used lineages → {args.out_used}")
    if no_score14_species:
        pd.DataFrame({"species": sorted(no_score14_species)}).to_csv(args.out_no_score14, index=False)
        print(f"🧭 Wrote species without score-1/4 protein → {args.out_no_score14}")
    if args.tree_check and not_in_tree:
        pd.DataFrame({"species": sorted(not_in_tree)}).to_csv(args.out_not_in_tree, index=False)
        print(f"🌱 Wrote matched-but-not-in-tree species → {args.out_not_in_tree}")

    # --- Build Newick from used lineages ---
    if used_rows:
        used_pairs = [(sp, lineage) for (sp, lineage, *_rest) in used_rows]
        newick = build_tree_from_used(used_pairs, verbose=args.verbose)
        Path(args.out_newick).write_text(newick, encoding="utf-8")
        print(f"🌳 Wrote Newick → {args.out_newick}")
    else:
        print("❌ No matched organisms; Newick not created.")

    print("✅ Done.")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import re
import pandas as pd
from pathlib import Path
from collections import Counter

# ---- Config ----
USED_CSV = "used_lineages_from_score14.csv"  # columns: species,lineage,source, n_entries_matched_1, n_entries_matched_4
OUT_NEWICK = "taxonomy_tree_sanitized.nwk"
DROP_SUPERROOTS = {"cellular organisms"}
USE_SPECIES_AS_LEAVES = True
BRANCH_LENGTH = 1

# ---- Helpers ----
def split_lineage(s: str, sep_guess=None):
    if not isinstance(s, str):
        return []
    s = s.strip()
    if not s:
        return []
    if sep_guess is None:
        sep = "," if s.count(",") >= s.count(";") else ";"
    else:
        sep = sep_guess
    parts = [p.strip() for p in s.split(sep)]
    return [p for p in parts if p]

def strip_rank(tag: str):
    if not isinstance(tag, str):
        return ""
    # drop ' (rank)' tails
    m = re.match(r"^(.*?)(?:\s*\([^)]*\))\s*$", tag)
    name = m.group(1) if m else tag
    return re.sub(r"\s+", " ", name).strip()

def clean_taxon_name(name: str):
    name = strip_rank(name)
    # normalize spaces & quotes
    name = name.replace("\xa0", " ").strip().strip('"').strip("'")
    return name

def lineage_to_nodes(lineage: str):
    parts = [clean_taxon_name(p) for p in split_lineage(lineage)]
    parts = [p for p in parts if p and p not in DROP_SUPERROOTS]
    # remove accidental immediate repeats
    dedup = []
    for p in parts:
        if not dedup or dedup[-1] != p:
            dedup.append(p)
    return dedup

class TreeNode:
    __slots__ = ("name", "children", "is_leaf_tag")
    def __init__(self, name=None):
        self.name = name
        self.children = {}   # name -> TreeNode
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

# Strict label escape: keep letters/digits/underscore/period/pipe/hyphen
SAFE_CHARS_RE = re.compile(r"[^A-Za-z0-9_.|\-]+")
def escape_newick_label(label: str) -> str:
    if label is None:
        return ""
    label = label.replace(" ", "_")
    # Replace any problematic char (including , : ; ( ) [ ] { } / etc.)
    label = SAFE_CHARS_RE.sub("_", label)
    # Newick disallows empty quoted names in many parsers; ensure not empty
    return label if label else "unnamed"

def to_newick(node: TreeNode, is_root=False):
    """Render Newick; omit branch length on root."""
    if not node.children:
        label = escape_newick_label(node.name) if node.name is not None else ""
        return f"{label}:{BRANCH_LENGTH}"
    subs = ",".join(to_newick(child, is_root=False) for child in node.children.values())
    label = "" if node.name is None else escape_newick_label(node.name)
    if is_root:
        # many parsers prefer no branch length on the root
        return f"({subs}){label}"
    else:
        return f"({subs}){label}:{BRANCH_LENGTH}"

def rebuild_newick_from_used(used_csv: str, out_nwk: str):
    path = Path(used_csv)
    if not path.exists():
        raise FileNotFoundError(f"Used-lineages CSV not found: {used_csv}")

    df = pd.read_csv(path, dtype=str)
    # Accept either "species, lineage" or "organism, lineage" column names
    cols = {c.lower(): c for c in df.columns}
    sp_col = cols.get("species") or cols.get("organism")
    lin_col = cols.get("lineage")
    if not sp_col or not lin_col:
        raise ValueError(f"{used_csv} must have columns including 'species' (or 'organism') and 'lineage'")

    pairs = []
    for _, row in df.iterrows():
        sp = (row.get(sp_col) or "").strip()
        lin = (row.get(lin_col) or "").strip()
        if sp and lin:
            pairs.append((sp, lin))
    if not pairs:
        raise ValueError("No (species, lineage) pairs found in used CSV")

    # Build tree
    root = TreeNode(name=None)
    dup = Counter()
    attached = 0
    for sp, lin in pairs:
        path = lineage_to_nodes(lin)
        if not path:
            continue
        if USE_SPECIES_AS_LEAVES:
            leaf = sp
            dup[leaf] += 1
            if dup[leaf] > 1:
                leaf = f"{leaf}__{dup[sp]}"
            root.add_path(path, leaf_label=leaf)
        else:
            root.add_path(path)
        attached += 1

    newick = to_newick(root, is_root=True)
    if not newick.endswith(";"):
        newick += ";"
    Path(out_nwk).write_text(newick, encoding="utf-8")
    print(f"🌳 Rebuilt sanitized Newick → {out_nwk} (attached {attached} leaves)")

if __name__ == "__main__":
    rebuild_newick_from_used(USED_CSV, OUT_NEWICK)

# #!/usr/bin/env python3
# import pandas as pd
# from pathlib import Path

# # === CONFIG ===
# PRED_DIR = Path("prediction_files")
# OUT_ORG_MAX = "organism_max_annotation.csv"
# OUT_SUMMARY = "organism_max_annotation_summary.csv"

# # --- 1) Load all prediction files (Organism, Entry, Annotation) ---
# csv_files = sorted(PRED_DIR.glob("*.csv"))
# if not csv_files:
#     raise FileNotFoundError(f"No CSV files found in {PRED_DIR}")

# frames = []
# for p in csv_files:
#     try:
#         df = pd.read_csv(p, dtype=str, on_bad_lines="skip")
#     except Exception as e:
#         print(f"⚠️ Could not read {p.name}: {e}")
#         continue
#     if not {"Organism", "Annotation"}.issubset(df.columns):
#         print(f"⚠️ Missing Organism/Annotation in {p.name}, skipping")
#         continue
#     df = df[["Organism", "Annotation"]].dropna()
#     frames.append(df)

# if not frames:
#     raise ValueError("No valid prediction files with Organism + Annotation columns")

# all_df = pd.concat(frames, ignore_index=True)
# all_df["Annotation"] = pd.to_numeric(all_df["Annotation"], errors="coerce").fillna(0).astype(int)
# print(f"📂 Loaded {len(all_df):,} rows across {len(csv_files)} files")

# # --- 2) Group by Organism and get max annotation ---
# org_max = all_df.groupby("Organism")["Annotation"].max().reset_index()
# print(f"🧬 Computed max annotations for {len(org_max):,} organisms")

# # --- 3) Count how many organisms max at each annotation value ---
# summary = org_max["Annotation"].value_counts().sort_index()
# print("\n=== Organism max annotation counts ===")
# print(summary)

# # --- 4) Save outputs ---
# org_max.to_csv(OUT_ORG_MAX, index=False)
# summary.to_csv(OUT_SUMMARY, header=["count"])
# print(f"\n📝 Wrote detailed per-organism max annotations → {OUT_ORG_MAX}")
# print(f"📝 Wrote summary counts → {OUT_SUMMARY}")

