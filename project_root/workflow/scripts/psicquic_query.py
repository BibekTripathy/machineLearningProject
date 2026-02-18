import pandas as pd
import os
import re
import subprocess
import mygene

# 1. Load config and input
psicquic_path = os.path.join(os.getcwd(), snakemake.config["psicquic"]["local_db"]["path"])
min_score = snakemake.config["psicquic"].get("min_score", 0.0)

seeds = pd.read_csv(snakemake.input["seeds"], sep="\t")
genes = sorted(seeds["gene_symbol"].unique().tolist())

print(f"Querying LOCAL IntAct/PSICQUIC database for {len(genes)} genes and their neighbors...")

# 2. Map Symbols to Uniprot IDs for a broader search
print("Mapping symbols to Uniprot IDs...")
mg = mygene.MyGeneInfo()
mapping = mg.querymany(genes, scopes="symbol", fields="uniprot", species=9606)

uniprot_ids = set()
for item in mapping:
    if "uniprot" in item:
        u = item["uniprot"]
        if "Swiss-Prot" in u:
            if isinstance(u["Swiss-Prot"], list): uniprot_ids.update(u["Swiss-Prot"])
            else: uniprot_ids.add(u["Swiss-Prot"])
        elif "TrEMBL" in u:
            if isinstance(u["TrEMBL"], list): uniprot_ids.update(u["TrEMBL"])
            else: uniprot_ids.add(u["TrEMBL"])

print(f"Total search identifiers: {len(genes)} symbols + {len(uniprot_ids)} Uniprot IDs")

# 3. Extract relevant rows efficiently using GREP
tmp_output = "results/networks/tmp_grep_psicquic.tsv"
os.makedirs("results/networks", exist_ok=True)

# Broad pattern: Match any of the symbols or Uniprot IDs
# We escape potential special characters and use word boundaries if possible, 
# but for MITAB, simple OR is usually best.
search_terms = set(genes) | uniprot_ids
pattern_file = "results/networks/grep_patterns.txt"
with open(pattern_file, "w") as f:
    for term in search_terms:
        f.write(f"{term}\n")

# Use fgrep (fixed strings) with -f (file) for massive speed boost with many patterns
grep_cmd = f"fgrep -i -f '{pattern_file}' '{psicquic_path}' > '{tmp_output}' || true"

print(f"Running broad fgrep extraction...")
subprocess.run(grep_cmd, shell=True, check=True)

if not os.path.exists(tmp_output) or os.path.getsize(tmp_output) == 0:
    print("Warning: No interactions found even with broad search.")
    direct_df = pd.DataFrame(columns=["nodeA", "nodeB", "weight", "source", "version"])
    expanded_df = pd.DataFrame(columns=["nodeA", "nodeB", "weight", "source", "version"])
else:
    # 4. Process the results
    print("Processing extracted interactions...")
    # Read the grepped rows
    df = pd.read_csv(tmp_output, sep="\t", header=None, low_memory=False)
    
    def parse_score(s):
        try:
            match = re.search(r":(\d+\.\d+|\d+)", str(s))
            return float(match.group(1)) if match else 0.5
        except: return 0.5
    
    def extract_symbol(s):
        # IntAct MITAB columns 0, 1 (IDs) and 4, 5 (Aliases)
        # We try to find a gene symbol first, then fall back to ID
        s_str = str(s)
        # Look for "(gene name)" or "(display_long)" or similar
        match = re.search(r":(.*?)\(", s_str)
        if match: return match.group(1)
        return s_str.split(":")[-1]

    # Map to our 40 genes specifically if they appear in aliases
    # This helps normalize the "neighbors" back to symbols where possible
    
    output_df = pd.DataFrame()
    # In MITAB 2.7: 0=ID A, 1=ID B, 4=Alias A, 5=Alias B, 14=Scores
    output_df["nodeA"] = df[4].apply(extract_symbol)
    output_df["nodeB"] = df[5].apply(extract_symbol)
    output_df["weight"] = df[14].apply(parse_score)
    output_df["source"] = "IntAct_Local_Broad"
    output_df["version"] = "MITAB_2.7"

    # Filter by score
    output_df = output_df[output_df["weight"] >= min_score]
    output_df = output_df[output_df["nodeA"] != output_df["nodeB"]]
    output_df = output_df.drop_duplicates()
    
    # --- Split into Direct and Expanded ---
    gene_set = set(genes)
    
    # Confirm nodes actually match our seeds or are real neighbors
    # (Since fgrep might match random strings in other columns)
    # We re-check the IDs as well (cols 0 and 1)
    df_ids_a = df[0].apply(lambda x: str(x).split(":")[-1])
    df_ids_b = df[1].apply(lambda x: str(x).split(":")[-1])
    
    # A row is valid if either Symbol OR Uniprot ID of either node is in our seed set
    valid_mask = (output_df["nodeA"].isin(gene_set) | df_ids_a.isin(uniprot_ids) |
                  output_df["nodeB"].isin(gene_set) | df_ids_b.isin(uniprot_ids))
    
    expanded_df = output_df[valid_mask].copy()
    
    # Direct: Both sides must be in the seed set (Symbol or Uniprot)
    direct_mask = ( (output_df["nodeA"].isin(gene_set) | df_ids_a.isin(uniprot_ids)) & 
                    (output_df["nodeB"].isin(gene_set) | df_ids_b.isin(uniprot_ids)) )
    direct_df = output_df[direct_mask].copy()

    # Cleanup
    if os.path.exists(tmp_output): os.remove(tmp_output)
    if os.path.exists(pattern_file): os.remove(pattern_file)

# 5. Save both
os.makedirs(os.path.dirname(snakemake.output["net"]), exist_ok=True)
direct_df.to_csv(snakemake.output["net"], sep="\t", index=False)
print(f"Direct PPI: {len(direct_df)} edges")
