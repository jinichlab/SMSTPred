#!/usr/bin/env python3
"""
Extract InterPro IDs for each Edison trend cluster and produce a frequency table.

For each cluster: lists unique InterPro IDs in order of frequency,
along with the total protein count per cluster.
"""

import pandas as pd
import numpy as np
import requests
import time
import json
from pathlib import Path
from collections import Counter
from tqdm import tqdm

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR = Path(__file__).parent.resolve()
PRED_DIR = SCRIPT_DIR.parent.parent / "classifier" / "reference_files" / "prediction_files"

# Edison trend clusters
EDISON_CLUSTERS = {
    "Vertebrate-Protein": [126, 186, 180, 185, 121],
    "Fungi-Other": [42, 167, 168, 14, 85, 169, 166],
    "Actino-Other": [26, 15, 5, 35, 198],
    "DNA-RM": [154, 139, 105, 107, 112],
    "Vertebrate-DNA": [48],
    "Parasite-RNA": [157, 205, 204],
    "Plant": [12, 18, 37],
}

# Flat set of all Edison cluster IDs
ALL_EDISON_CLUSTER_IDS = set()
for cids in EDISON_CLUSTERS.values():
    ALL_EDISON_CLUSTER_IDS.update(cids)

# ============================================================================
# Main
# ============================================================================

def main():
    # Load cluster assignments
    cluster_df = pd.read_csv(SCRIPT_DIR / "cluster_assignments.csv")
    print(f"Loaded {len(cluster_df)} entries with cluster assignments")

    # Filter to only Edison trend cluster entries
    edison_mask = cluster_df['cluster'].isin(ALL_EDISON_CLUSTER_IDS)
    edison_entries = set(cluster_df.loc[edison_mask, 'Entry'].values)
    print(f"Edison trend cluster entries: {len(edison_entries)}")

    # Build Entry -> InterPro mapping from prediction files
    print(f"\nScanning {len(list(PRED_DIR.glob('*_predictions.csv')))} prediction files for InterPro IDs...")
    entry_interpro = {}

    for pred_file in tqdm(list(PRED_DIR.glob("*_predictions.csv")), desc="Reading prediction files"):
        try:
            df = pd.read_csv(pred_file, usecols=["Entry", "InterPro"], low_memory=False)
            for _, row in df.iterrows():
                entry = row["Entry"]
                if entry in edison_entries and entry not in entry_interpro:
                    interpro = row.get("InterPro", "")
                    if pd.notna(interpro) and str(interpro).strip():
                        entry_interpro[entry] = str(interpro).strip()
        except (ValueError, KeyError):
            # File might not have InterPro column
            continue

    print(f"Found InterPro annotations for {len(entry_interpro)} / {len(edison_entries)} entries")

    # Add InterPro to cluster dataframe
    cluster_df['InterPro'] = cluster_df['Entry'].map(entry_interpro)

    # Build frequency tables per cluster
    all_rows = []

    for trend_name, cluster_ids in EDISON_CLUSTERS.items():
        for cid in cluster_ids:
            cdf = cluster_df[cluster_df['cluster'] == cid]
            total_proteins = len(cdf)

            if total_proteins == 0:
                continue

            # Count InterPro IDs
            interpro_counter = Counter()
            proteins_with_interpro = 0

            for _, row in cdf.iterrows():
                ipr = row.get('InterPro', '')
                if pd.notna(ipr) and str(ipr).strip():
                    proteins_with_interpro += 1
                    ids = [x.strip() for x in str(ipr).split(';') if x.strip()]
                    for ipr_id in ids:
                        interpro_counter[ipr_id] += 1

            # Sort by frequency
            for ipr_id, count in interpro_counter.most_common():
                all_rows.append({
                    'Trend': trend_name,
                    'Cluster': cid,
                    'Total_Proteins': total_proteins,
                    'Proteins_With_InterPro': proteins_with_interpro,
                    'InterPro_ID': ipr_id,
                    'Count': count,
                    'Frequency': count / total_proteins,
                })

    result_df = pd.DataFrame(all_rows)

    # ========================================================================
    # Fetch InterPro descriptions from EBI API
    # ========================================================================
    unique_ipr_ids = result_df['InterPro_ID'].unique()
    print(f"\nFetching descriptions for {len(unique_ipr_ids)} unique InterPro IDs from EBI API...")

    ipr_descriptions = {}
    ipr_types = {}
    cache_path = SCRIPT_DIR / "interpro_description_cache.json"

    # Load cache if it exists
    if cache_path.exists():
        with open(cache_path) as f:
            cache = json.load(f)
        print(f"  Loaded {len(cache)} cached descriptions")
    else:
        cache = {}

    for ipr_id in tqdm(unique_ipr_ids, desc="Fetching InterPro descriptions"):
        if ipr_id in cache:
            ipr_descriptions[ipr_id] = cache[ipr_id].get('name', 'N/A')
            ipr_types[ipr_id] = cache[ipr_id].get('type', 'N/A')
            continue

        url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/{ipr_id}"
        try:
            resp = requests.get(url, headers={"Accept": "application/json"}, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                name = data.get("metadata", {}).get("name", {})
                if isinstance(name, dict):
                    name = name.get("name", "N/A")
                entry_type = data.get("metadata", {}).get("type", "N/A")
                ipr_descriptions[ipr_id] = name
                ipr_types[ipr_id] = entry_type
                cache[ipr_id] = {'name': name, 'type': entry_type}
            elif resp.status_code == 404:
                ipr_descriptions[ipr_id] = "Not found"
                ipr_types[ipr_id] = "N/A"
                cache[ipr_id] = {'name': 'Not found', 'type': 'N/A'}
            else:
                ipr_descriptions[ipr_id] = f"HTTP {resp.status_code}"
                ipr_types[ipr_id] = "N/A"
            time.sleep(0.15)  # Rate limit: ~6 requests/sec
        except Exception as e:
            ipr_descriptions[ipr_id] = f"Error: {e}"
            ipr_types[ipr_id] = "N/A"

    # Save cache
    with open(cache_path, 'w') as f:
        json.dump(cache, f, indent=2)
    print(f"  Saved description cache ({len(cache)} entries)")

    # Add descriptions to result DataFrame
    result_df['InterPro_Name'] = result_df['InterPro_ID'].map(ipr_descriptions)
    result_df['InterPro_Type'] = result_df['InterPro_ID'].map(ipr_types)

    # Save full table
    output_path = SCRIPT_DIR / "interpro_by_cluster.csv"
    result_df.to_csv(output_path, index=False)
    print(f"\nSaved full table to: {output_path}")

    # Print summary
    print("\n" + "=" * 100)
    print("INTERPRO IDs BY EDISON TREND CLUSTER")
    print("=" * 100)

    for trend_name, cluster_ids in EDISON_CLUSTERS.items():
        print(f"\n{'─' * 100}")
        print(f"TREND: {trend_name}")
        print(f"{'─' * 100}")

        for cid in cluster_ids:
            tdf = result_df[result_df['Cluster'] == cid]
            if len(tdf) == 0:
                continue

            total = tdf['Total_Proteins'].iloc[0]
            with_ipr = tdf['Proteins_With_InterPro'].iloc[0]
            print(f"\n  Cluster {cid} (n={total}, {with_ipr} with InterPro):")
            print(f"  {'InterPro ID':<12} {'Count':>5} {'Freq':>7}  {'Type':<12} {'Name'}")
            print(f"  {'─'*85}")

            for _, row in tdf.head(15).iterrows():
                name = str(row.get('InterPro_Name', 'N/A'))[:50]
                itype = str(row.get('InterPro_Type', 'N/A'))
                print(f"  {row['InterPro_ID']:<12} {row['Count']:>5} {row['Frequency']:>6.1%}  {itype:<12} {name}")

            if len(tdf) > 15:
                print(f"  ... and {len(tdf) - 15} more InterPro IDs")

    print(f"\n\nFull results saved to: {output_path}")


if __name__ == "__main__":
    main()
