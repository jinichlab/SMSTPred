import pandas as pd
import re
import requests
import time
import io
from pathlib import Path

# ========= CONFIGURATION =========
INPUT_FILE = "../reference_files/reference_proteomes.csv"  # from your filtered .csv
OUTPUT_DIR = Path("../reference_files/reference_proteomes")
OUTPUT_DIR.mkdir(exist_ok=True)
API_SLEEP = 1.5  # seconds between API calls to avoid rate limits

# ========= FILTERING =========
def get_busco_completeness(busco_string):
    if pd.isna(busco_string):
        return None
    match = re.search(r"C:(\d+\.\d+)%", busco_string)
    return float(match.group(1)) if match else None

print("🔍 Reading input file...")
df = pd.read_csv(INPUT_FILE)
df['BUSCO_completeness'] = df['BUSCO'].apply(get_busco_completeness)

# Apply filters
df = df[~df['Taxonomic lineage'].str.contains('Viruses', na=False)].copy()
df = df[df['BUSCO_completeness'].fillna(0) >= 95].copy()
df = df.loc[df.groupby('Organism Id')['Protein count'].idxmax()]
print(f"✅ {len(df)} proteomes after filtering.")

# ========= UNIPROT API QUERY =========
fields = "accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Csequence%2Cannotation_score%2Cxref_interpro%2Cgo_id%2Cec%2Crhea%2Clit_pubmed_id%2Cxref_pfam%2Ccc_function"
base_url = "https://rest.uniprot.org/uniprotkb/stream"

for idx, row in df.iterrows():
    proteome_id = row['Proteome Id']
    label = row['Organism'].replace(" ", "_").replace("/", "_")
    out_csv = OUTPUT_DIR / f"{label}.csv"

    print(f"\n⏳ Fetching proteome {proteome_id} → {out_csv.name}")
    url = f"{base_url}?fields={fields}&format=tsv&query=(proteome:{proteome_id})"

    try:
        response = requests.get(url)
        response.raise_for_status()

        proteins_df = pd.read_csv(io.StringIO(response.text), sep='\t')
        proteins_df['Proteome Id'] = proteome_id

        original_count = len(proteins_df)
        proteins_df = proteins_df[proteins_df["Length"] <= 1000]
        print(f"🧪 Filtered {original_count - len(proteins_df)} proteins >1000aa")

        # Write selected columns
        subset = proteins_df[[
            'Entry', 'Protein names', 'Organism', 'Length', 'Sequence', 'InterPro', 'Pfam'
        ]]
        subset.to_csv(out_csv, index=False)

        print(f"✅ Saved {len(subset)} proteins to {out_csv.name}")
    
    except requests.exceptions.RequestException as e:
        print(f"❌ Error fetching {proteome_id}: {e}")

    time.sleep(API_SLEEP)  # avoid hammering UniProt