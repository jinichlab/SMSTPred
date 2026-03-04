#!/usr/bin/env python3
"""
Visualize UMAP clusters colored by cluster ID, with Edison trend clusters highlighted.
Reviewed (SwissProt) entries are shown with distinct markers for easy lookup.

This script creates visualizations of the clustered SAM MTase UMAP projection,
highlighting clusters that match the Edison AI discoveries.

Requirements:
    pip install datamapplot pandas numpy matplotlib seaborn requests
"""

import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import colorsys
import gzip
import io

# Try to import datamapplot for interactive visualization
try:
    import datamapplot
    HAS_DATAMAPPLOT = True
except ImportError:
    HAS_DATAMAPPLOT = False
    print("Warning: datamapplot not installed. Only static plots will be generated.")

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR = Path(__file__).parent.resolve()
OUTPUT_DIR = SCRIPT_DIR
PREDICTOR_DIR = SCRIPT_DIR.parent.parent / "substrate_type_predictor_testing"
PRED_DIR = PREDICTOR_DIR / "reference_files" / "prediction_files"

# Key clusters identified from Edison trend analysis
EDISON_CLUSTERS = {
    # Trend 1: Vertebrate + Protein
    "Vertebrate-Protein": {
        "clusters": [126, 186, 180, 185, 121],
        "color": "#2ecc71",  # Green
        "description": "Vertebrate protein methylation"
    },
    # Trend 2: Fungi + Other (secondary metabolism)
    "Fungi-Other": {
        "clusters": [42, 167, 168, 14, 85, 169, 166],
        "color": "#e67e22",  # Orange
        "description": "Fungi secondary metabolism"
    },
    # Trend 3: Actinomycetota + Other (BGC)
    "Actino-Other": {
        "clusters": [26, 15, 5, 35, 198],
        "color": "#9b59b6",  # Purple
        "description": "Actinomycetota BGC methylation"
    },
    # Trend 4: DNA clusters (RM systems)
    "DNA-RM": {
        "clusters": [154, 139, 105, 107, 112],
        "color": "#3498db",  # Blue
        "description": "DNA methylation (RM systems)"
    },
    # Trend 5: Vertebrate DNA (epigenetics)
    "Vertebrate-DNA": {
        "clusters": [48],
        "color": "#1abc9c",  # Teal
        "description": "Vertebrate DNA methylation (DNMTs)"
    },
    # Trend 6: Parasite + RNA
    "Parasite-RNA": {
        "clusters": [157, 205, 204],
        "color": "#e74c3c",  # Red
        "description": "Parasite RNA modification"
    },
    # Trend 7: Plant clusters
    "Plant": {
        "clusters": [12, 18, 37],
        "color": "#27ae60",  # Dark green
        "description": "Plant-specific MTases"
    },
}

# Create a flat mapping of cluster -> trend name
CLUSTER_TO_TREND = {}
for trend_name, info in EDISON_CLUSTERS.items():
    for cid in info["clusters"]:
        CLUSTER_TO_TREND[cid] = trend_name


def generate_cluster_colors(n_clusters):
    """Generate distinct colors for clusters using HSV color space."""
    colors = {}
    for i in range(n_clusters):
        hue = i / n_clusters
        saturation = 0.6 + (i % 3) * 0.15
        value = 0.7 + (i % 2) * 0.2
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        colors[i] = '#%02x%02x%02x' % (int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
    return colors


def load_reviewed_entries_from_file():
    """
    Load set of reviewed (SwissProt) entries from local Excel file.

    The file contains all Swiss-Prot reviewed SAM MTases downloaded from UniProt.
    Any entry present in this file is considered reviewed.
    """
    reviewed_file = SCRIPT_DIR / "uniprotkb_IPR029063_AND_reviewed_true_2026_02_05.xlsx.gz"

    if not reviewed_file.exists():
        raise FileNotFoundError(
            f"Reviewed entries file not found: {reviewed_file}\n"
            "Please download Swiss-Prot reviewed SAM MTases from UniProt."
        )

    print(f"Loading reviewed entries from: {reviewed_file.name}")

    # Decompress gzipped Excel file and read with openpyxl engine
    with gzip.open(reviewed_file, 'rb') as f:
        reviewed_df = pd.read_excel(io.BytesIO(f.read()), engine='openpyxl')

    # Get set of reviewed entry IDs
    reviewed_entries = set(reviewed_df['Entry'].unique())
    print(f"  Found {len(reviewed_entries):,} reviewed (SwissProt) entries")

    return reviewed_entries


def load_cluster_data():
    """Load cluster assignments and original UMAP data."""
    cluster_path = SCRIPT_DIR / "cluster_assignments.csv"

    if not cluster_path.exists():
        raise FileNotFoundError(f"Cluster assignments not found at {cluster_path}")

    df = pd.read_csv(cluster_path)
    print(f"Loaded {len(df)} entries with cluster assignments")
    print(f"Clusters found: {df['cluster'].nunique()} (including noise)")

    # Load reviewed entries from local file
    reviewed_entries = load_reviewed_entries_from_file()

    # Add reviewed column - True if entry is in the reviewed set
    df['Reviewed'] = df['Entry'].isin(reviewed_entries)
    n_reviewed = df['Reviewed'].sum()
    print(f"Reviewed (SwissProt) entries in dataset: {n_reviewed:,} ({n_reviewed/len(df)*100:.1f}%)")

    return df


def create_hover_text(df):
    """Create hover text for each point."""
    hover_texts = []

    for _, row in df.iterrows():
        cluster = int(row['cluster'])
        trend = CLUSTER_TO_TREND.get(cluster, "")
        trend_str = f"\n*** {trend} ***" if trend else ""
        reviewed_str = "✓ SwissProt (Reviewed)" if row.get('Reviewed', False) else "TrEMBL (Unreviewed)"

        text = (
            f"Entry: {row['Entry']}\n"
            f"Status: {reviewed_str}\n"
            f"Organism: {row['Organism']}\n"
            f"Category: {row['Organism_category']}\n"
            f"Substrate: {row['Substrate']}\n"
            f"Cluster: {cluster}{trend_str}"
        )
        hover_texts.append(text)

    return np.array(hover_texts)


def get_cluster_color(cluster_id, cluster_colors):
    """Get color for a cluster, prioritizing Edison trend colors."""
    if cluster_id == -1:
        return "#cccccc"

    trend = CLUSTER_TO_TREND.get(cluster_id)
    if trend:
        return EDISON_CLUSTERS[trend]["color"]

    return cluster_colors.get(cluster_id, "#888888")


def create_static_plot(df, output_path):
    """Create static matplotlib visualization colored by cluster, with reviewed entries highlighted."""
    fig, axes = plt.subplots(1, 2, figsize=(24, 10))

    clusters = sorted([c for c in df['cluster'].unique() if c != -1])
    cluster_colors = generate_cluster_colors(len(clusters))

    # -------------------------------------------------------------------------
    # Left plot: All clusters, reviewed entries with stars
    # -------------------------------------------------------------------------
    ax1 = axes[0]

    # Plot noise first
    noise_mask = df['cluster'] == -1
    if noise_mask.sum() > 0:
        ax1.scatter(
            df.loc[noise_mask, 'UMAP_1'],
            df.loc[noise_mask, 'UMAP_2'],
            c='#dddddd',
            alpha=0.3,
            s=3,
            marker='.'
        )

    # Plot each cluster - unreviewed first, then reviewed on top
    for cluster_id in clusters:
        mask = df['cluster'] == cluster_id
        color = get_cluster_color(cluster_id, cluster_colors)
        is_edison = cluster_id in CLUSTER_TO_TREND

        # Unreviewed entries (circles)
        unrev_mask = mask & ~df['Reviewed']
        if unrev_mask.sum() > 0:
            ax1.scatter(
                df.loc[unrev_mask, 'UMAP_1'],
                df.loc[unrev_mask, 'UMAP_2'],
                c=color,
                alpha=0.5,
                s=5,
                marker='o',
                edgecolors='none',
            )

        # Reviewed entries (stars, larger, with edge)
        rev_mask = mask & df['Reviewed']
        if rev_mask.sum() > 0:
            ax1.scatter(
                df.loc[rev_mask, 'UMAP_1'],
                df.loc[rev_mask, 'UMAP_2'],
                c=color,
                alpha=0.9,
                s=60 if is_edison else 40,
                marker='*',
                edgecolors='black',
                linewidths=0.5,
            )

    ax1.set_xlabel("UMAP 1", fontsize=12)
    ax1.set_ylabel("UMAP 2", fontsize=12)
    ax1.set_title("All Clusters (★ = SwissProt Reviewed)", fontsize=14, fontweight="bold")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Add legend for marker types
    legend_elements = [
        plt.scatter([], [], c='gray', s=40, marker='*', edgecolors='black', linewidths=0.5, label='Reviewed (SwissProt)'),
        plt.scatter([], [], c='gray', s=20, marker='o', edgecolors='none', label='Unreviewed (TrEMBL)'),
    ]
    ax1.legend(handles=legend_elements, loc='upper right', fontsize=10)

    # -------------------------------------------------------------------------
    # Right plot: Only Edison trend clusters, reviewed highlighted
    # -------------------------------------------------------------------------
    ax2 = axes[1]

    # Plot all points as gray background
    ax2.scatter(
        df['UMAP_1'],
        df['UMAP_2'],
        c='#eeeeee',
        alpha=0.2,
        s=2,
    )

    legend_handles = []

    for trend_name, info in EDISON_CLUSTERS.items():
        trend_mask = df['cluster'].isin(info['clusters'])
        if trend_mask.sum() > 0:
            # Unreviewed
            unrev_mask = trend_mask & ~df['Reviewed']
            if unrev_mask.sum() > 0:
                ax2.scatter(
                    df.loc[unrev_mask, 'UMAP_1'],
                    df.loc[unrev_mask, 'UMAP_2'],
                    c=info['color'],
                    alpha=0.6,
                    s=10,
                    marker='o',
                    edgecolors='none',
                )

            # Reviewed (stars)
            rev_mask = trend_mask & df['Reviewed']
            if rev_mask.sum() > 0:
                ax2.scatter(
                    df.loc[rev_mask, 'UMAP_1'],
                    df.loc[rev_mask, 'UMAP_2'],
                    c=info['color'],
                    alpha=0.95,
                    s=100,
                    marker='*',
                    edgecolors='black',
                    linewidths=0.5,
                )

            n_reviewed = rev_mask.sum()
            n_total = trend_mask.sum()
            legend_handles.append(
                mpatches.Patch(
                    color=info['color'],
                    label=f"{trend_name}: {info['description']} (★{n_reviewed}/{n_total})"
                )
            )

            # Add cluster number labels
            for cid in info['clusters']:
                cluster_mask = df['cluster'] == cid
                if cluster_mask.sum() > 0:
                    cx = df.loc[cluster_mask, 'UMAP_1'].mean()
                    cy = df.loc[cluster_mask, 'UMAP_2'].mean()
                    ax2.annotate(
                        str(cid),
                        (cx, cy),
                        fontsize=8,
                        fontweight='bold',
                        ha='center',
                        va='center',
                        color='white',
                        bbox=dict(boxstyle='circle,pad=0.3', facecolor=info['color'], edgecolor='black', linewidth=0.5)
                    )

    ax2.set_xlabel("UMAP 1", fontsize=12)
    ax2.set_ylabel("UMAP 2", fontsize=12)
    ax2.set_title("Edison Trend Clusters (★ = SwissProt Reviewed)", fontsize=14, fontweight="bold")
    ax2.legend(handles=legend_handles, loc='upper left', fontsize=9, framealpha=0.9)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"Saved static plot to: {output_path}")


def create_cluster_detail_plots(df, output_dir):
    """Create individual plots for each Edison trend category, highlighting reviewed entries."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    for trend_name, info in EDISON_CLUSTERS.items():
        fig, ax = plt.subplots(figsize=(12, 10))

        # Gray background
        ax.scatter(df['UMAP_1'], df['UMAP_2'], c='#eeeeee', alpha=0.2, s=2)

        reviewed_entries = []

        for cid in info['clusters']:
            cluster_mask = df['cluster'] == cid
            if cluster_mask.sum() > 0:
                cluster_df = df[cluster_mask]

                # Unreviewed (circles)
                unrev = cluster_df[~cluster_df['Reviewed']]
                if len(unrev) > 0:
                    ax.scatter(
                        unrev['UMAP_1'], unrev['UMAP_2'],
                        c=info['color'], alpha=0.6, s=15,
                        marker='o', edgecolors='none',
                        label=f"Cluster {cid} unreviewed (n={len(unrev)})"
                    )

                # Reviewed (stars)
                rev = cluster_df[cluster_df['Reviewed']]
                if len(rev) > 0:
                    ax.scatter(
                        rev['UMAP_1'], rev['UMAP_2'],
                        c=info['color'], alpha=0.95, s=150,
                        marker='*', edgecolors='black', linewidths=0.8,
                        label=f"Cluster {cid} reviewed (n={len(rev)})"
                    )
                    reviewed_entries.extend(rev['Entry'].tolist())

                # Cluster label
                cx, cy = cluster_df['UMAP_1'].mean(), cluster_df['UMAP_2'].mean()
                ax.annotate(
                    str(cid), (cx, cy),
                    fontsize=12, fontweight='bold', ha='center', va='center',
                    color='white',
                    bbox=dict(boxstyle='circle,pad=0.4', facecolor=info['color'], edgecolor='black', linewidth=1)
                )

        ax.set_xlabel("UMAP 1", fontsize=12)
        ax.set_ylabel("UMAP 2", fontsize=12)
        ax.set_title(f"{trend_name}: {info['description']}\n(★ = {len(reviewed_entries)} SwissProt entries)", fontsize=14, fontweight="bold")
        ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plot_path = output_dir / f"cluster_{trend_name.lower().replace('-', '_')}.png"
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
        plt.close()

        print(f"Saved {trend_name} detail plot to: {plot_path}")

        # Save list of reviewed entries for this trend
        if reviewed_entries:
            entries_path = output_dir / f"reviewed_entries_{trend_name.lower().replace('-', '_')}.txt"
            with open(entries_path, 'w') as f:
                f.write(f"# Reviewed (SwissProt) entries for {trend_name}\n")
                f.write(f"# {info['description']}\n")
                f.write(f"# Clusters: {info['clusters']}\n\n")
                for entry in reviewed_entries:
                    f.write(f"{entry}\n")
            print(f"  Saved {len(reviewed_entries)} reviewed entries to: {entries_path}")


def create_interactive_plot(df, output_path):
    """Create interactive datamapplot visualization with reviewed status."""
    if not HAS_DATAMAPPLOT:
        print("Skipping interactive plot (datamapplot not installed)")
        return

    clusters = sorted([c for c in df['cluster'].unique() if c != -1])
    cluster_colors = generate_cluster_colors(len(clusters))

    # Color array - reviewed entries are gold for clear visibility
    base_colors = df['cluster'].apply(lambda c: get_cluster_color(c, cluster_colors)).values
    color_array = np.where(df['Reviewed'], '#FFD700', base_colors)

    # Labels include reviewed status
    def get_label(row):
        cid = int(row['cluster'])
        reviewed = "★ REVIEWED" if row.get('Reviewed', False) else ""
        if cid == -1:
            return f"Noise {reviewed}"
        trend = CLUSTER_TO_TREND.get(cid, "")
        if trend:
            return f"Cluster {cid} ({trend}) {reviewed}"
        return f"Cluster {cid} {reviewed}"

    labels = df.apply(get_label, axis=1).values

    # Marker sizes - much larger for reviewed so they stand out
    marker_sizes = np.where(df['Reviewed'], 25, 5)

    data_map = df[["UMAP_1", "UMAP_2"]].to_numpy().astype("float32")
    hover_data = create_hover_text(df)

    custom_css = """
    .row { display: flex; align-items: center; }
    .box { height:12px; width:12px; border-radius:3px; margin-right:8px; }
    #legend {
        position: absolute; bottom: 0; right: 0; margin: 16px; padding: 14px;
        border-radius: 12px; z-index: 2; background: #ffffffdd;
        font-family: Arial, sans-serif; font-size: 10pt;
        box-shadow: 2px 3px 10px #aaaaaa44; max-height: 400px; overflow-y: auto;
    }
    #title-container { max-width: 75%; }
    """

    custom_html = '<div id="legend"><b>Edison Trend Clusters</b><br><br>\n'
    custom_html += f'<div class="row"><div class="box" style="background-color:#FFD700; border: 1px solid #B8960F;"></div><b>Gold = SwissProt Reviewed</b></div><br>\n'
    for trend_name, info in EDISON_CLUSTERS.items():
        custom_html += f'<div class="row"><div class="box" style="background-color:{info["color"]};"></div>{trend_name}: {", ".join(map(str, info["clusters"]))}</div>\n'
    custom_html += "</div>"

    print("Creating interactive plot...")
    plot = datamapplot.create_interactive_plot(
        data_map,
        labels,
        font_family="Arial",
        title="SAM MTase UMAP - Cluster Analysis",
        sub_title=f"Gold = SwissProt Reviewed | {len(df):,} MTases | {df['Reviewed'].sum():,} reviewed",
        hover_text=hover_data,
        enable_search=True,
        marker_color_array=color_array,
        custom_css=custom_css,
        custom_html=custom_html,
        marker_size_array=marker_sizes,
    )

    plot.save(str(output_path))
    print(f"Saved interactive plot to: {output_path}")


def print_cluster_summary(df):
    """Print summary of Edison trend clusters with reviewed counts."""
    print("\n" + "="*70)
    print("EDISON TREND CLUSTER SUMMARY")
    print("="*70)

    for trend_name, info in EDISON_CLUSTERS.items():
        trend_mask = df['cluster'].isin(info['clusters'])
        n_total = trend_mask.sum()
        n_reviewed = (trend_mask & df['Reviewed']).sum()

        print(f"\n{trend_name}: {info['description']}")
        print(f"  Total MTases: {n_total:,} (★ {n_reviewed} reviewed)")
        print(f"  Clusters: {info['clusters']}")

        if n_total > 0:
            trend_df = df[trend_mask]
            print(f"  Substrate breakdown:")
            for sub in ['RNA', 'DNA', 'Protein', 'Other']:
                sub_pct = (trend_df['Substrate'] == sub).sum() / n_total * 100
                print(f"    {sub}: {sub_pct:.1f}%")

            # List reviewed entries
            rev_entries = trend_df[trend_df['Reviewed']]['Entry'].tolist()
            if rev_entries:
                print(f"  Reviewed entries ({len(rev_entries)}):")
                for entry in rev_entries[:10]:
                    print(f"    - {entry}")
                if len(rev_entries) > 10:
                    print(f"    ... and {len(rev_entries)-10} more")


def main():
    print("="*70)
    print("Cluster Visualization for Edison Trend Analysis")
    print("(with SwissProt Reviewed entries highlighted)")
    print("="*70)

    # Load data
    print("\n[1/4] Loading cluster data and reviewed status...")
    df = load_cluster_data()

    # Print summary
    print_cluster_summary(df)

    # Create output paths
    date_str = datetime.now().strftime("%Y-%m-%d")

    # Static plot
    print("\n[2/4] Creating static overview plot...")
    static_path = OUTPUT_DIR / f"clusters_overview_{date_str}.png"
    create_static_plot(df, static_path)

    # Detail plots
    print("\n[3/4] Creating trend detail plots...")
    detail_dir = OUTPUT_DIR / "cluster_details"
    create_cluster_detail_plots(df, detail_dir)

    # Interactive plot
    print("\n[4/4] Creating interactive plot...")
    interactive_path = OUTPUT_DIR / f"clusters_interactive_{date_str}.html"
    create_interactive_plot(df, interactive_path)

    # Summary
    print("\n" + "="*70)
    print("COMPLETE!")
    print("="*70)
    print(f"\nOutput files:")
    print(f"  - {static_path.name} (overview)")
    print(f"  - cluster_details/ (individual trend plots + reviewed entry lists)")
    if HAS_DATAMAPPLOT:
        print(f"  - {interactive_path.name} (interactive HTML)")
    print("\n★ = SwissProt Reviewed entries (good for literature lookup)")
    print("Open the HTML file in a browser for interactive exploration!")


if __name__ == "__main__":
    main()
