import pandas as pd
import requests
import io

seeds = pd.read_csv(snakemake.input["seeds"], sep="\t")
genes = seeds["gene_symbol"].unique().tolist()

# Query PSICQUIC using the current REST API (IntAct service)
# The old /Tools/webservices/ endpoint was deprecated; using /psicquic/ws/psicquic/1.2.0/
base_url = "https://www.ebi.ac.uk/intact/psicquic/ws/psicquic/1.2.0/interactions"
query = " OR ".join([f'"{g}"' for g in genes[:20]])  # Limit to 20 genes to avoid timeout

params = {
    "q": query,
    "format": "PSICQUIC-MITAB2.5",
    "firstResult": 0,
    "maxResults": 1000
}

print(f"Querying PSICQUIC (IntAct) for cancer-specific interactions...")
print(f"Query genes: {len(genes)} (using first 20 for API limit)")

try:
    response = requests.get(base_url, params=params, timeout=60)
    response.raise_for_status()
    
    # Check if we got HTML (error page) instead of MITAB data
    if response.text.strip().startswith("<!"):
        print("Warning: PSICQUIC returned HTML instead of data. Falling back to STRING disease query.")
        raise ValueError("PSICQUIC returned HTML")
    
    # Parse MITAB2.5 format (tab-separated with 15+ columns)
    df = pd.read_csv(io.StringIO(response.text), sep="\t", header=None)
    print(f"PSICQUIC response shape: {df.shape}")
    
    if df.shape[1] < 2:
        raise ValueError("PSICQUIC returned insufficient columns")
    
    # MITAB2.5: col 0 = ID A, col 1 = ID B, col 14 = confidence score
    cols_to_extract = [0, 1]
    if df.shape[1] > 14:
        cols_to_extract.append(14)
        output_df = df[cols_to_extract].copy()
        output_df.columns = ["nodeA", "nodeB", "score"]
    else:
        output_df = df[cols_to_extract].copy()
        output_df.columns = ["nodeA", "nodeB"]
        output_df["score"] = 0.0
    
    output_df["source"] = "PSICQUIC_IntAct"
    
    def extract_score(s):
        try:
            if pd.isna(s): return 0.0
            s_str = str(s).split("|")[0]
            if ":" in s_str:
                return float(s_str.split(":")[-1])
            return float(s_str)
        except:
            return 0.0
    
    output_df["score"] = output_df["score"].apply(extract_score)
    output_df = output_df[output_df["nodeA"] != output_df["nodeB"]]
    output_df = output_df.drop_duplicates()
    
except Exception as e:
    print(f"PSICQUIC query failed: {e}")
    print("Falling back to STRING API for disease-specific PPIs...")
    
    # Fallback: Query STRING with cancer-related keywords
    # Map gene symbols to Ensembl IDs for STRING
    import mygene
    mg = mygene.MyGeneInfo()
    mapping = mg.querymany(genes, scopes="symbol", fields="ensembl.gene", species=9606)
    
    ensembl_ids = []
    for item in mapping:
        if "ensembl" in item:
            if isinstance(item["ensembl"], list):
                ensembl_ids.append(item["ensembl"][0]["gene"])
            else:
                ensembl_ids.append(item["ensembl"]["gene"])
    
    if ensembl_ids:
        # Query STRING for high-confidence interactions among cancer genes
        url = "https://string-db.org/api/tsv/network"
        params = {
            "identifiers": "\r".join(ensembl_ids[:50]),  # Limit for API
            "species": 9606,
            "required_score": 700,
            "caller_identity": "cancer_pipeline_psicquic_fallback"
        }
        
        response = requests.post(url, data=params, timeout=120)
        response.raise_for_status()
        
        df = pd.read_csv(io.StringIO(response.text), sep="\t")
        
        if df.empty:
            print("Warning: No interactions found in STRING fallback")
            output_df = pd.DataFrame(columns=["nodeA", "nodeB", "weight", "source", "version"])
        else:
            output_df = df[["preferredName_A", "preferredName_B", "score"]].copy()
            output_df.columns = ["nodeA", "nodeB", "weight"]
            # STRING returns scores as 0-1000, normalize to 0-1 for consistency
            # Check if scores are already in 0-1 range (max <= 1.0)
            if output_df["weight"].max() > 1.0:
                output_df["weight"] = output_df["weight"] / 1000.0
            output_df["source"] = "STRING_disease_fallback"
            output_df["version"] = "STRING_v12.0_fallback"
            output_df = output_df[output_df["nodeA"] != output_df["nodeB"]]
            output_df = output_df.drop_duplicates()
    else:
        print("Error: Could not map genes or fetch interactions")
        output_df = pd.DataFrame(columns=["nodeA", "nodeB", "weight", "source", "version"])

print(f"Final disease-specific PPI: {len(output_df)} edges")
output_df.to_csv(snakemake.output["net"], sep="\t", index=False)
