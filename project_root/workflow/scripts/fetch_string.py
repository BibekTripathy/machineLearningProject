import pandas as pd
import os

# 1. Load config and input
links_path = os.path.join(os.getcwd(), snakemake.config["string"]["local_db"]["links_path"])
info_path = os.path.join(os.getcwd(), snakemake.config["string"]["local_db"]["info_path"])
min_score = snakemake.config["string"]["min_combined_score"] * 1000 # STRING file uses 0-1000

seeds = pd.read_csv(snakemake.input[0], sep="\t")
symbols = seeds["gene_symbol"].unique().tolist()

print(f"Loading local STRING data for {len(symbols)} seed genes...")

# 2. Map Symbols to STRING IDs (using protein.info)
info_df = pd.read_csv(info_path, sep="\t")
# The info file has columns: #string_protein_id, preferred_name, protein_size, annotation
# We map preferred_name (Symbol) to string_protein_id (ENSP)
mapping_df = info_df[info_df["preferred_name"].isin(symbols)]
symbol_to_id = dict(zip(mapping_df["preferred_name"], mapping_df["#string_protein_id"]))
id_to_symbol = dict(zip(mapping_df["#string_protein_id"], mapping_df["preferred_name"]))

string_ids = list(symbol_to_id.values())

if not string_ids:
    raise ValueError("None of the seed genes could be mapped to STRING IDs in the local info file.")

print(f"Mapped {len(string_ids)} symbols to STRING IDs.")

# 3. Filter STRING links
# The links file can be very large, so we use chunking or efficient filtering if possible
# Since we only care about interactions where BOTH nodes are in our set (or at least one?)
# Usually, for a background network, we might want all interactions involving our seeds.
# For this pipeline, we'll fetch interactions where AT LEAST ONE node is a seed, 
# and the other is also in our seed list (or expanded list? No, Phase 2A is the backbone).
# Actually, the original API call fetched the network for those identifiers.
# Let's filter for interactions where both protein1 and protein2 are in our seed list
# to replicate the "network" behavior for the initial backbone.

print(f"Filtering local links file: {links_path}")
# Note: For very large files, pd.read_csv(chunksize=...) is better.
# For 9606, it's ~600MB, which fits in memory on most systems.
links_df = pd.read_csv(links_path, sep=" ")

# Filter by score
links_df = links_df[links_df["combined_score"] >= min_score]

# Filter for our seed genes
# Option A: Both must be in seeds (strictly local network)
# Option B: One must be in seeds (seeds + neighbors)
# The original script used the STRING API 'network' endpoint which returns interactions among the input genes.
mask = links_df["protein1"].isin(string_ids) & links_df["protein2"].isin(string_ids)
filtered_links = links_df[mask].copy()

# 4. Format and Save
version = snakemake.config["string"]["version_tag"]

if filtered_links.empty:
    print("Warning: No interactions found for the given genes at this threshold in local DB.")
    output_df = pd.DataFrame(columns=["nodeA", "nodeB", "weight", "source", "evidence", "version"])
else:
    # Map back to symbols
    filtered_links["nodeA"] = filtered_links["protein1"].map(id_to_symbol)
    filtered_links["nodeB"] = filtered_links["protein2"].map(id_to_symbol)
    filtered_links["weight"] = filtered_links["combined_score"] / 1000.0
    
    output_df = filtered_links[["nodeA", "nodeB", "weight"]].copy()
    output_df["source"] = "STRING_Local"
    output_df["evidence"] = "combined_score"
    output_df["version"] = version

    # Remove duplicates and self-loops
    output_df = output_df[output_df["nodeA"] != output_df["nodeB"]]
    output_df = output_df.drop_duplicates()

os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)
output_df.to_csv(snakemake.output[0], sep="\t", index=False)
print(f"Local network backbone saved with {len(output_df)} edges.")
