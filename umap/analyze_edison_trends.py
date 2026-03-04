#!/usr/bin/env python3
"""
Analyze UMAP clusters to find patterns matching Edison AI discoveries:
1. Vertebrate protein-binding enrichment
2. Free-living bacteria "Other" enrichment
3. Helicobacter DNA specialization
4. Parasite RNA specialization (Toxoplasma)
5. Fungi/Actinomycetota secondary metabolism

Run from: /home/nroethler/Code/UMAP/sam_umap_output/
"""

import pandas as pd
import numpy as np
from sklearn.cluster import HDBSCAN
from collections import Counter

# Load data
df = pd.read_csv('sam_mtase_umap_2026-02-04.csv')

print(f"Total entries: {len(df)}")
print(f"\nSubstrate distribution:")
print(df['Substrate'].value_counts())

# Define organism categories
vertebrate_keywords = [
    'homo sapiens', 'human', 'mus musculus', 'mouse', 'rattus', 'rat',
    'danio rerio', 'zebrafish', 'gallus gallus', 'chicken',
    'bos taurus', 'cattle', 'sus scrofa', 'pig', 'canis lupus', 'dog',
    'xenopus', 'frog', 'cyprinus carpio', 'carp', 'oryzias', 'medaka',
    'heterocephalus glaber', 'naked mole rat', 'capra hircus', 'goat',
    'ovis aries', 'sheep', 'equus', 'horse', 'felis catus', 'cat',
    'macaca', 'pan troglodytes', 'chimpanzee', 'gorilla', 'cirrhinus',
    'salmo', 'oncorhynchus', 'takifugu', 'tetraodon', 'gasterosteus',
    'anolis', 'chelonia', 'alligator', 'crocodylus', 'python', 'boa'
]

bacteria_keywords = [
    'bacterium', 'bacillus', 'coccus', 'monas', 'bacter', 'vibrio',
    'clostridium', 'escherichia', 'salmonella', 'pseudomonas',
    'lactobacillus', 'staphylococcus', 'enterococcus', 'klebsiella',
    'acinetobacter', 'burkholderia', 'rhizobium', 'agrobacterium',
    'caulobacter', 'sphingomonas', 'xanthomonas', 'erwinia',
    'geobacter', 'shewanella', 'desulfovibrio', 'methan', 'nostoc',
    'synechocystis', 'prochlorococcus', 'cyanobact'
]

actinomycetota_keywords = [
    'streptomyces', 'mycobacterium', 'amycolatopsis', 'nocardia',
    'corynebacterium', 'rhodococcus', 'frankia', 'micromonospora',
    'actinomyces', 'bifidobacterium', 'propionibacterium', 'geodermatophilus'
]

fungi_keywords = [
    'fungus', 'yeast', 'aspergillus', 'fusarium', 'candida',
    'saccharomyces', 'penicillium', 'cryptococcus', 'neurospora',
    'trichoderma', 'botrytis', 'magnaporthe', 'ustilago', 'puccinia',
    'schizosaccharomyces', 'yarrowia', 'pichia', 'kluyveromyces',
    'colletotrichum', 'verticillium', 'metarhizium', 'beauveria',
    'glarea', 'zasmidium', 'aureobasidium', 'rhodotorula', 'moniliophthora',
    'phytophthora', 'pythium', 'rhizophagus', 'glomus', 'coccomyxa'
]

plant_keywords = [
    'arabidopsis', 'oryza', 'rice', 'zea mays', 'maize', 'triticum', 'wheat',
    'populus', 'lotus', 'nelumbo', 'prunus', 'pomegranate', 'punica',
    'cotton', 'gossypium', 'solanum', 'nicotiana', 'vitis', 'citrus',
    'cucumis', 'citrullus', 'brassica', 'glycine', 'soybean', 'medicago',
    'phaseolus', 'ziziphus', 'malus', 'pyrus', 'fragaria', 'rosa'
]

parasite_keywords = [
    'toxoplasma', 'plasmodium', 'trypanosoma', 'leishmania', 'giardia',
    'cryptosporidium', 'entamoeba', 'trichomonas', 'babesia', 'theileria',
    'eimeria', 'neospora', 'cyclospora'
]

helicobacter_keywords = ['helicobacter', 'campylobacter']

def classify_organism(org_name):
    """Classify organism into broad categories"""
    if pd.isna(org_name):
        return 'Unknown'
    org_lower = org_name.lower()

    # Check specific groups first (order matters)
    for h in helicobacter_keywords:
        if h in org_lower:
            return 'Helicobacter'

    for p in parasite_keywords:
        if p in org_lower:
            return 'Parasite'

    for v in vertebrate_keywords:
        if v in org_lower:
            return 'Vertebrate'

    for a in actinomycetota_keywords:
        if a in org_lower:
            return 'Actinomycetota'

    for f in fungi_keywords:
        if f in org_lower:
            return 'Fungi'

    for p in plant_keywords:
        if p in org_lower:
            return 'Plant'

    for b in bacteria_keywords:
        if b in org_lower:
            return 'Bacteria'

    return 'Other_organism'

df['Organism_category'] = df['Organism'].apply(classify_organism)

print(f"\nOrganism category distribution:")
print(df['Organism_category'].value_counts())

# Cluster using HDBSCAN
print("\n" + "="*80)
print("CLUSTERING WITH HDBSCAN")
print("="*80)

coords = df[['UMAP_1', 'UMAP_2']].values
clusterer = HDBSCAN(min_cluster_size=50, min_samples=10)
df['cluster'] = clusterer.fit_predict(coords)

n_clusters = len(set(df['cluster'])) - (1 if -1 in df['cluster'].values else 0)
print(f"\nFound {n_clusters} clusters (excluding noise)")
print(f"Noise points: {(df['cluster'] == -1).sum()}")

# Analyze each cluster
def analyze_cluster(cluster_df, cluster_id):
    """Analyze a single cluster for composition"""
    n = len(cluster_df)
    substrate_dist = cluster_df['Substrate'].value_counts(normalize=True) * 100
    org_dist = cluster_df['Organism_category'].value_counts(normalize=True) * 100
    centroid = (cluster_df['UMAP_1'].mean(), cluster_df['UMAP_2'].mean())
    return {
        'n': n,
        'substrate_dist': substrate_dist,
        'org_dist': org_dist,
        'centroid': centroid,
        'top_organisms': cluster_df['Organism'].value_counts().head(5)
    }

cluster_analyses = {}
for cluster_id in sorted(df['cluster'].unique()):
    if cluster_id == -1:
        continue
    cluster_df = df[df['cluster'] == cluster_id]
    cluster_analyses[cluster_id] = analyze_cluster(cluster_df, cluster_id)

# Find interesting clusters matching Edison trends
print("\n" + "="*80)
print("EDISON TREND ANALYSIS")
print("="*80)

print("\n" + "-"*80)
print("TREND 1: Vertebrate Protein-binding enrichment")
print("-"*80)

found_vert_protein = False
for cid, analysis in cluster_analyses.items():
    vert_pct = analysis['org_dist'].get('Vertebrate', 0)
    protein_pct = analysis['substrate_dist'].get('Protein', 0)
    if vert_pct > 25 and protein_pct > 40:
        found_vert_protein = True
        print(f"\n*** CLUSTER {cid} - Vertebrate Protein cluster ***")
        print(f"Size: {analysis['n']} MTases")
        print(f"Centroid: ({analysis['centroid'][0]:.2f}, {analysis['centroid'][1]:.2f})")
        print(f"Vertebrate: {vert_pct:.1f}%")
        print(f"Protein-binding: {protein_pct:.1f}%")
        print(f"Top organisms:\n{analysis['top_organisms']}")
if not found_vert_protein:
    print("No clusters found with >25% vertebrate AND >40% protein-binding")
    # Show best candidates
    best = max(cluster_analyses.items(),
               key=lambda x: x[1]['org_dist'].get('Vertebrate', 0) * x[1]['substrate_dist'].get('Protein', 0))
    print(f"Best candidate: Cluster {best[0]} with {best[1]['org_dist'].get('Vertebrate', 0):.1f}% vertebrate, {best[1]['substrate_dist'].get('Protein', 0):.1f}% protein")

print("\n" + "-"*80)
print("TREND 2: Fungi/Actinomycetota 'Other' (secondary metabolism)")
print("-"*80)

found_fungi_other = False
for cid, analysis in cluster_analyses.items():
    fungi_pct = analysis['org_dist'].get('Fungi', 0)
    actino_pct = analysis['org_dist'].get('Actinomycetota', 0)
    other_pct = analysis['substrate_dist'].get('Other', 0)
    if (fungi_pct > 30 or actino_pct > 15) and other_pct > 50:
        found_fungi_other = True
        print(f"\n*** CLUSTER {cid} - Fungi/Actinomycetota secondary metabolism ***")
        print(f"Size: {analysis['n']} MTases")
        print(f"Centroid: ({analysis['centroid'][0]:.2f}, {analysis['centroid'][1]:.2f})")
        print(f"Fungi: {fungi_pct:.1f}%, Actinomycetota: {actino_pct:.1f}%")
        print(f"Other-binding: {other_pct:.1f}%")
        print(f"Top organisms:\n{analysis['top_organisms']}")
if not found_fungi_other:
    print("No clusters found with high fungi/actinomycetota AND >50% Other-binding")

print("\n" + "-"*80)
print("TREND 3: DNA-enriched clusters (Helicobacter specialization)")
print("-"*80)

found_dna = False
for cid, analysis in cluster_analyses.items():
    dna_pct = analysis['substrate_dist'].get('DNA', 0)
    helico_pct = analysis['org_dist'].get('Helicobacter', 0)
    if dna_pct > 60:
        found_dna = True
        print(f"\n*** CLUSTER {cid} - DNA-enriched cluster ***")
        print(f"Size: {analysis['n']} MTases")
        print(f"Centroid: ({analysis['centroid'][0]:.2f}, {analysis['centroid'][1]:.2f})")
        print(f"DNA-binding: {dna_pct:.1f}%")
        print(f"Helicobacter: {helico_pct:.1f}%")
        print(f"Organism breakdown:\n{analysis['org_dist'].head()}")
        print(f"Top organisms:\n{analysis['top_organisms']}")
if not found_dna:
    print("No clusters found with >60% DNA-binding")

print("\n" + "-"*80)
print("TREND 4: RNA-enriched clusters (check for parasite enrichment)")
print("-"*80)

found_rna = False
for cid, analysis in cluster_analyses.items():
    rna_pct = analysis['substrate_dist'].get('RNA', 0)
    parasite_pct = analysis['org_dist'].get('Parasite', 0)
    if rna_pct > 70:
        found_rna = True
        print(f"\n*** CLUSTER {cid} - RNA-enriched cluster ***")
        print(f"Size: {analysis['n']} MTases")
        print(f"Centroid: ({analysis['centroid'][0]:.2f}, {analysis['centroid'][1]:.2f})")
        print(f"RNA-binding: {rna_pct:.1f}%")
        print(f"Parasite: {parasite_pct:.1f}%")
        print(f"Organism breakdown:\n{analysis['org_dist'].head()}")
        print(f"Top organisms:\n{analysis['top_organisms']}")
if not found_rna:
    print("No clusters found with >70% RNA-binding")

# Summary of all clusters
print("\n" + "="*80)
print("SUMMARY OF ALL CLUSTERS")
print("="*80)

for cid in sorted(cluster_analyses.keys()):
    analysis = cluster_analyses[cid]
    dominant_substrate = analysis['substrate_dist'].idxmax()
    dominant_org = analysis['org_dist'].idxmax()
    print(f"Cluster {cid:2d}: n={analysis['n']:4d}, "
          f"Substrate: {dominant_substrate:8s} ({analysis['substrate_dist'].max():.0f}%), "
          f"Organism: {dominant_org:15s} ({analysis['org_dist'].max():.0f}%), "
          f"Centroid: ({analysis['centroid'][0]:6.2f}, {analysis['centroid'][1]:6.2f})")

# Save cluster assignments
output_path = 'cluster_assignments.csv'
output_df = df[['Entry', 'Organism', 'Substrate', 'UMAP_1', 'UMAP_2', 'cluster', 'Organism_category']]
output_df.to_csv(output_path, index=False)
print(f"\nCluster assignments saved to {output_path}")

# Also save a detailed report for each interesting cluster
print("\n" + "="*80)
print("DETAILED CLUSTER REPORTS")
print("="*80)

for cid, analysis in cluster_analyses.items():
    print(f"\n--- Cluster {cid} ---")
    print(f"Size: {analysis['n']}")
    print(f"Centroid: ({analysis['centroid'][0]:.2f}, {analysis['centroid'][1]:.2f})")
    print(f"\nSubstrate distribution:")
    for sub, pct in analysis['substrate_dist'].items():
        print(f"  {sub}: {pct:.1f}%")
    print(f"\nOrganism category distribution:")
    for org, pct in analysis['org_dist'].items():
        print(f"  {org}: {pct:.1f}%")
    print(f"\nTop 5 organisms:")
    for org, count in analysis['top_organisms'].items():
        print(f"  {org}: {count}")
