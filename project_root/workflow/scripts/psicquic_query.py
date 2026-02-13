import pandas as pd
import requests
import io

seeds = pd.read_csv(snakemake.input["seeds"], sep="\t")
genes = seeds["gene_symbol"].unique().tolist()

# 2. Query PSICQUIC (Example: IntAct service)
# We use MIQL query language as per Phase 3A [cite: 45]
base_url = "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/"
query = " OR ".join([f"identifier:{g}" for g in genes])
full_url = f"{base_url}{query}"

print(f"Querying PSICQUIC for cancer-specific interactions...")
response = requests.get(full_url)

# 3. Format into an edgelist [cite: 53]
# PSICQUIC returns Tab25 format; we extract the gene names
df = pd.read_csv(io.StringIO(response.text), sep="\t", header=None)
# Simplification: extracting column 1 and 2 (Gene IDs)
output_df = df[[0, 1]].copy()
output_df.columns = ["nodeA", "nodeB"]
output_df["source"] = "PSICQUIC_IntAct"

output_df.to_csv(snakemake.output["net"], sep="\t", index=False)
