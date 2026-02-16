import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import re

# Load both sets of modules
bg_mod = pd.read_csv(snakemake.input["bg_mod"], sep="\t")
dis_mod = pd.read_csv(snakemake.input["dis_mod"], sep="\t")

# Extract numeric module_id for comparison (remove 'leiden_' prefix)
def extract_module_num(mod_id):
    if isinstance(mod_id, str):
        match = re.search(r'leiden_(\d+)', mod_id)
        if match:
            return int(match.group(1))
    return mod_id

bg_mod['module_num'] = bg_mod['module_id'].apply(extract_module_num)
dis_mod['module_num'] = dis_mod['module_id'].apply(extract_module_num)

# We only compare genes that are present in both clustering results
common_genes = list(set(bg_mod["gene_symbol"]) & set(dis_mod["gene_symbol"]))

if not common_genes:
    ari = 0.0
    nmi = 0.0
    jaccard = 0.0
else:
    # Filter and sort to ensure labels match
    bg_sub = bg_mod[bg_mod["gene_symbol"].isin(common_genes)].sort_values("gene_symbol")
    dis_sub = dis_mod[dis_mod["gene_symbol"].isin(common_genes)].sort_values("gene_symbol")

    ari = adjusted_rand_score(bg_sub["module_num"], dis_sub["module_num"])
    nmi = normalized_mutual_info_score(bg_sub["module_num"], dis_sub["module_num"])
    
    # Calculate Jaccard similarity of module assignments
    bg_sets = bg_sub.groupby("module_num")["gene_symbol"].apply(set).to_dict()
    dis_sets = dis_sub.groupby("module_num")["gene_symbol"].apply(set).to_dict()
    
    # Average Jaccard across modules
    jaccards = []
    for mod_id, genes in bg_sets.items():
        for dis_mod_id, dis_genes in dis_sets.items():
            if len(genes.union(dis_genes)) > 0:
                jaccard = len(genes.intersection(dis_genes)) / len(genes.union(dis_genes))
                jaccards.append(jaccard)
    
    jaccard = sum(jaccards) / len(jaccards) if jaccards else 0.0

# Save metrics with additional ML-relevant information
metrics = pd.DataFrame({
    "Metric": ["Adjusted Rand Index (ARI)", "Normalized Mutual Information (NMI)", "Average Jaccard Similarity"],
    "Value": [ari, nmi, jaccard],
    "Description": [
        "Agreement between clustering assignments (0=random, 1=perfect)",
        "Mutual information between clusterings (0=none, 1=perfect)",
        "Average overlap of genes across modules"
    ]
})
metrics.to_csv(snakemake.output["metrics"], sep="\t", index=False)
