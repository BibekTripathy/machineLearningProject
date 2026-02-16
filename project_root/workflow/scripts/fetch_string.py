import pandas as pd
import mygene
import requests
import io
import os

# 1. Initialize MyGene and load input
mg = mygene.MyGeneInfo()

# Load seeds (your 40 genes) from the path defined in the Snakefile
seeds = pd.read_csv(snakemake.input[0], sep="\t")
symbols = seeds["gene_symbol"].unique().tolist()

# 2. Map Symbols to Ensembl IDs (Standardizing to GRCh38)
print(f"Mapping {len(symbols)} symbols to Ensembl IDs...")
mapping = mg.querymany(symbols, scopes="symbol", fields="ensembl.gene", species=9606)

# Extract Ensembl IDs, filtering out any genes that didn't map
ensembl_ids = []
for item in mapping:
    if "ensembl" in item:
        # Some genes might have multiple Ensembl IDs; we take the first one
        if isinstance(item["ensembl"], list):
            ensembl_ids.append(item["ensembl"][0]["gene"])
        else:
            ensembl_ids.append(item["ensembl"]["gene"])

if not ensembl_ids:
    raise ValueError("No genes were successfully mapped to Ensembl IDs.")

# 3. Fetch STRING interactions using the mapped IDs [cite: 303, 305]
# Using 'network' endpoint for functional/physical associations [cite: 13, 165]
url = "https://string-db.org/api/tsv/network"
params = {
    "identifiers": "\r".join(ensembl_ids),
    "species": 9606,  # Homo sapiens [cite: 229]
    "required_score": 700,  # 0.7 threshold scaled to 1000 [cite: 230]
    "caller_identity": "cancer_stratification_pipeline",
}

print(f"Fetching STRING interactions for {len(ensembl_ids)} mapped IDs...")
response = requests.post(url, data=params)
response.raise_for_status()

# 4. Process and Format [cite: 310]
# We keep Node A, Node B, and the combined score [cite: 101, 312]
df = pd.read_csv(io.StringIO(response.text), sep="\t")

# Get version from config
version = snakemake.config.get("string", {}).get("version_tag", "v12.0")

# Handle empty responses
if df.empty:
    print("Warning: No interactions found for the given genes at this threshold.")
    output_df = pd.DataFrame(columns=["nodeA", "nodeB", "weight", "source", "evidence", "version"])
else:
    output_df = df[["preferredName_A", "preferredName_B", "score"]].copy()
    output_df.columns = ["nodeA", "nodeB", "weight"]

    # STRING TSV format returns scores as 0-1000, normalize to 0-1 for consistency
    # Check if scores are already in 0-1 range (max <= 1.0)
    if output_df["weight"].max() > 1.0:
        output_df["weight"] = output_df["weight"] / 1000.0

    # Add metadata columns per framework schema
    output_df["source"] = "STRING"
    output_df["evidence"] = "combined_score"
    output_df["version"] = version

    # Remove duplicates and self-loops [cite: 311, 312]
    output_df = output_df[output_df["nodeA"] != output_df["nodeB"]]
    output_df = output_df.drop_duplicates()

# 5. Save Output to the path defined in Snakefile
os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)
output_df.to_csv(snakemake.output[0], sep="\t", index=False)
print(f"Network backbone saved with {len(output_df)} edges.")
