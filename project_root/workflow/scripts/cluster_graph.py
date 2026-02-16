import pandas as pd
import igraph as ig
import leidenalg
import networkx as nx

# Determine disease type from output path
output_path = snakemake.output[0]
if "disease_modules" in output_path:
    disease = "Cancer_Disease_PPI"
else:
    disease = "Cancer_Background"

# Load the network
net = pd.read_csv(snakemake.input["net"], sep="\t")

# Determine which column has the weight (could be 'score' or 'weight')
weight_col = "weight" if "weight" in net.columns else "score"

if "expanded" in snakemake.input.keys():
    # Load high-ranking genes from RWR
    expanded = pd.read_csv(snakemake.input["expanded"], sep="\t")
    genes = set(expanded["gene_symbol"])
    # Filter network to only include expanded genes
    sub_net = net[net["nodeA"].isin(genes) & net["nodeB"].isin(genes)].copy()
else:
    sub_net = net.copy()

# Convert to igraph for Leiden algorithm [cite: 321]
# Handle both 'weight' and 'score' column names
if weight_col in sub_net.columns:
    g = ig.Graph.TupleList(
        sub_net[["nodeA", "nodeB", weight_col]].itertuples(index=False),
        edge_attrs=[weight_col],
        weights=False  # weights are specified in edge_attrs
    )
else:
    g = ig.Graph.TupleList(
        sub_net[["nodeA", "nodeB"]].itertuples(index=False),
        weights=False
    )

# Run Leiden [cite: 42, 321]
partition = leidenalg.find_partition(
    g, leidenalg.ModularityVertexPartition, weights=weight_col
)

# Map results back to gene names
results = pd.DataFrame({
    "gene_symbol": g.vs["name"],
    "module_id": partition.membership,
    "method": "leiden",
    "disease": disease
})

# Calculate module quality metrics for ML feature selection
module_sizes = results.groupby("module_id").size().to_dict()
results["module_size"] = results["module_id"].map(module_sizes)

# Calculate module density using networkx
G_nx = nx.from_pandas_edgelist(sub_net, "nodeA", "nodeB", [weight_col])
module_density = {}
for mod_id in results["module_id"].unique():
    mod_genes = results[results["module_id"] == mod_id]["gene_symbol"].tolist()
    subgraph = G_nx.subgraph(mod_genes)
    n = len(subgraph.nodes())
    m = len(subgraph.edges())
    # Density = 2*m / (n*(n-1)) for undirected graphs
    if n > 1:
        density = 2 * m / (n * (n - 1))
    else:
        density = 0.0
    module_density[mod_id] = density

results["module_density"] = results["module_id"].map(module_density)

# Format module_id with method prefix for ML feature naming
results["module_id"] = results["module_id"].apply(lambda x: f"leiden_{x}")

results = results.sort_values(["module_id", "gene_symbol"]).reset_index(drop=True)
results.to_csv(snakemake.output[0], sep="\t", index=False)
