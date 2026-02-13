import pandas as pd
import networkx as nx
import numpy as np
from scipy import sparse

# 1. Load data
net = pd.read_csv(snakemake.input["net"], sep="\t")
seeds = pd.read_csv(snakemake.input["seeds"], sep="\t")

# 2. Build Graph
G = nx.from_pandas_edgelist(net, "nodeA", "nodeB", ["weight"])
nodes = sorted(G.nodes())
idx = {n: i for i, n in enumerate(nodes)}

# 3. Setup RWR [cite: 316, 317]
n = len(nodes)
p0 = np.zeros(n)
valid_seeds = [s for s in seeds["gene_symbol"] if s in idx]
if not valid_seeds:
    raise ValueError(
        "No seeds found in the network. Check gene symbols or network coverage."
    )
for s in valid_seeds:
    p0[idx[s]] = 1.0
p0 /= p0.sum()

# 4. Math: Iterative propagation
alpha = snakemake.config["propagation"]["rwr"]["restart_prob"]
A = nx.adjacency_matrix(G, nodelist=nodes)
# Normalize matrix columns to sum to 1
D_inv = sparse.diags(1.0 / np.array(A.sum(axis=0)).ravel())
T = A @ D_inv

p = p0.copy()
for _ in range(100):  # Iterations
    p = (1 - alpha) * (T @ p) + alpha * p0

# 5. Save Top K
top_idx = np.argsort(-p)[: snakemake.config["propagation"]["rwr"]["top_k"]]
results = pd.DataFrame(
    {"gene_symbol": [nodes[i] for i in top_idx], "score": p[top_idx]}
)
results.to_csv(snakemake.output[0], sep="\t", index=False)
