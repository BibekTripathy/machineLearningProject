import pandas as pd
import igraph as ig
import leidenalg

# Load the network and the high-ranking genes from RWR
net = pd.read_csv(snakemake.input["net"], sep="\t")
expanded = pd.read_csv(snakemake.input["expanded"], sep="\t")
genes = set(expanded["gene_symbol"])

# Filter network to only include expanded genes
sub_net = net[net["nodeA"].isin(genes) & net["nodeB"].isin(genes)]

# Convert to igraph for Leiden algorithm [cite: 321]
g = ig.Graph.TupleList(sub_net.itertuples(index=False), weights=True)

# Run Leiden [cite: 42, 321]
partition = leidenalg.find_partition(
    g, leidenalg.ModularityVertexPartition, weights="weight"
)

# Map results back to gene names
results = pd.DataFrame({"gene_symbol": g.vs["name"], "module_id": partition.membership})
results.to_csv(snakemake.output[0], sep="\t", index=False)
