# UMAP

UMAP dimensionality reduction and cluster analysis for SAM MTase embeddings.

## Scripts

| Script | Description |
|--------|-------------|
| `sam_mtase_substrate_umap.py` | Generates 2D UMAP projections of SAM MTases colored by predicted substrate type. Uses cosine metric, 15 neighbors, 0.1 min_dist. Outputs static PNG and interactive HTML (via datamapplot). |
| `analyze_edison_trends.py` | Analyzes trends in cluster composition across prediction iterations. |
| `quick_analysis.py` | Quick summary statistics and analysis of UMAP cluster results. |
| `visualize_clusters.py` | Cluster visualization with detailed annotation overlays. |
| `extract_interpro_by_cluster.py` | Extracts InterPro domain annotations for proteins within each UMAP cluster. |

## Usage

```bash
# Generate UMAP visualization
python sam_mtase_substrate_umap.py
```

## Output

- Static scatter plot (PNG) with substrate type coloring and legend
- Interactive HTML visualization with hover information and search
- CSV with UMAP coordinates and cluster assignments

## Environment

Uses the module-specific `environment.yml` in this directory (includes UMAP and datamapplot dependencies).
