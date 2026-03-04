#!/usr/bin/env python3
"""Quick analysis of cluster assignments to find Edison trends"""
import pandas as pd

df = pd.read_csv('cluster_assignments.csv')

# Get substrate info from original file
orig = pd.read_csv('sam_mtase_umap_2026-02-04.csv')
df['Substrate'] = orig['Substrate']

print("="*80)
print("CLUSTER ANALYSIS FOR EDISON TRENDS")
print("="*80)

# Analyze each cluster
results = []
for cid in df['cluster'].unique():
    if cid == -1:
        continue
    c = df[df['cluster'] == cid]
    n = len(c)

    # Substrate percentages
    sub_counts = c['Substrate'].value_counts(normalize=True) * 100
    rna_pct = sub_counts.get('RNA', 0)
    dna_pct = sub_counts.get('DNA', 0)
    protein_pct = sub_counts.get('Protein', 0)
    other_pct = sub_counts.get('Other', 0)

    # Organism category percentages
    org_counts = c['Organism_category'].value_counts(normalize=True) * 100
    vert_pct = org_counts.get('Vertebrate', 0)
    fungi_pct = org_counts.get('Fungi', 0)
    actino_pct = org_counts.get('Actinomycetota', 0)
    bacteria_pct = org_counts.get('Bacteria', 0)
    parasite_pct = org_counts.get('Parasite', 0)
    helico_pct = org_counts.get('Helicobacter', 0)
    plant_pct = org_counts.get('Plant', 0)

    centroid = (c['UMAP_1'].mean(), c['UMAP_2'].mean())

    results.append({
        'cluster': cid, 'n': n,
        'RNA': rna_pct, 'DNA': dna_pct, 'Protein': protein_pct, 'Other': other_pct,
        'Vertebrate': vert_pct, 'Fungi': fungi_pct, 'Actinomycetota': actino_pct,
        'Bacteria': bacteria_pct, 'Parasite': parasite_pct, 'Helicobacter': helico_pct,
        'Plant': plant_pct, 'centroid': centroid
    })

rdf = pd.DataFrame(results)

# TREND 1: Vertebrate + Protein
print("\n" + "="*80)
print("TREND 1: VERTEBRATE PROTEIN-BINDING CLUSTERS")
print("(Edison: Vertebrates expanded protein methylation)")
print("="*80)
vert_protein = rdf[(rdf['Vertebrate'] > 30) & (rdf['Protein'] > 40)].sort_values('Vertebrate', ascending=False)
if len(vert_protein) > 0:
    for _, row in vert_protein.head(10).iterrows():
        print(f"\nCluster {int(row['cluster'])}: n={int(row['n'])}, centroid=({row['centroid'][0]:.1f}, {row['centroid'][1]:.1f})")
        print(f"  Vertebrate: {row['Vertebrate']:.1f}%, Protein: {row['Protein']:.1f}%")
        c = df[df['cluster'] == row['cluster']]
        print(f"  Top organisms: {c['Organism'].value_counts().head(3).to_dict()}")
else:
    print("No clusters with >30% vertebrate AND >40% protein")
    # Show best candidates
    rdf['vert_prot_score'] = rdf['Vertebrate'] * rdf['Protein']
    best = rdf.nlargest(5, 'vert_prot_score')
    print("\nBest candidates:")
    for _, row in best.iterrows():
        print(f"  Cluster {int(row['cluster'])}: Vert={row['Vertebrate']:.1f}%, Prot={row['Protein']:.1f}%")

# TREND 2: Fungi/Actinomycetota + Other (secondary metabolism)
print("\n" + "="*80)
print("TREND 2: FUNGI/ACTINOMYCETOTA 'OTHER' CLUSTERS")
print("(Edison: BGC-linked secondary metabolism)")
print("="*80)
fungi_other = rdf[((rdf['Fungi'] > 25) | (rdf['Actinomycetota'] > 15)) & (rdf['Other'] > 50)].sort_values('Other', ascending=False)
if len(fungi_other) > 0:
    for _, row in fungi_other.head(10).iterrows():
        print(f"\nCluster {int(row['cluster'])}: n={int(row['n'])}, centroid=({row['centroid'][0]:.1f}, {row['centroid'][1]:.1f})")
        print(f"  Fungi: {row['Fungi']:.1f}%, Actino: {row['Actinomycetota']:.1f}%, Other: {row['Other']:.1f}%")
        c = df[df['cluster'] == row['cluster']]
        print(f"  Top organisms: {c['Organism'].value_counts().head(3).to_dict()}")
else:
    print("No clusters with high fungi/actino AND >50% other")

# TREND 3: DNA-enriched (Helicobacter/RM systems)
print("\n" + "="*80)
print("TREND 3: DNA-BINDING CLUSTERS")
print("(Edison: Helicobacter 52.5% DNA for restriction-modification defense)")
print("="*80)
dna_clusters = rdf[rdf['DNA'] > 50].sort_values('DNA', ascending=False)
if len(dna_clusters) > 0:
    for _, row in dna_clusters.head(10).iterrows():
        print(f"\nCluster {int(row['cluster'])}: n={int(row['n'])}, centroid=({row['centroid'][0]:.1f}, {row['centroid'][1]:.1f})")
        print(f"  DNA: {row['DNA']:.1f}%, Helicobacter: {row['Helicobacter']:.1f}%")
        c = df[df['cluster'] == row['cluster']]
        print(f"  Top organisms: {c['Organism'].value_counts().head(3).to_dict()}")
else:
    print("No clusters with >50% DNA")

# TREND 4: RNA + Parasite
print("\n" + "="*80)
print("TREND 4: RNA-BINDING CLUSTERS (check for parasite enrichment)")
print("(Edison: Toxoplasma 79.7% RNA for developmental regulation)")
print("="*80)
rna_clusters = rdf[rdf['RNA'] > 70].sort_values('RNA', ascending=False)
if len(rna_clusters) > 0:
    for _, row in rna_clusters.head(10).iterrows():
        print(f"\nCluster {int(row['cluster'])}: n={int(row['n'])}, centroid=({row['centroid'][0]:.1f}, {row['centroid'][1]:.1f})")
        print(f"  RNA: {row['RNA']:.1f}%, Parasite: {row['Parasite']:.1f}%")
        c = df[df['cluster'] == row['cluster']]
        print(f"  Top organisms: {c['Organism'].value_counts().head(3).to_dict()}")
else:
    print("No clusters with >70% RNA")

# TREND 5: Plant-specific clusters
print("\n" + "="*80)
print("TREND 5: PLANT-ENRICHED CLUSTERS")
print("="*80)
plant_clusters = rdf[rdf['Plant'] > 30].sort_values('Plant', ascending=False)
if len(plant_clusters) > 0:
    for _, row in plant_clusters.head(10).iterrows():
        print(f"\nCluster {int(row['cluster'])}: n={int(row['n'])}, centroid=({row['centroid'][0]:.1f}, {row['centroid'][1]:.1f})")
        print(f"  Plant: {row['Plant']:.1f}%, Dominant substrate: RNA={row['RNA']:.0f}%, DNA={row['DNA']:.0f}%, Prot={row['Protein']:.0f}%, Other={row['Other']:.0f}%")
        c = df[df['cluster'] == row['cluster']]
        print(f"  Top organisms: {c['Organism'].value_counts().head(3).to_dict()}")

# Save summary
rdf.to_csv('cluster_summary.csv', index=False)
print(f"\n\nFull cluster summary saved to cluster_summary.csv")
