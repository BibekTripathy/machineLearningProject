import pandas as pd
from gprofiler import (
    GProfiler,
)  # Ensure gprofiler-official is in your enrichment.yaml [cite: 60]

# Load modules and the background gene list
modules = pd.read_csv(snakemake.input["modules"], sep="\t")
background_net = pd.read_csv(snakemake.input["background"], sep="\t")
bg_genes = list(
    set(background_net["nodeA"]).union(set(background_net["nodeB"]))
)  # [cite: 59]

gp = GProfiler(return_dataframe=True)
all_enrichments = []

# Perform enrichment per module [cite: 289]
for mod_id in modules["module_id"].unique():
    query_genes = modules[modules["module_id"] == mod_id]["gene_symbol"].tolist()

    res = gp.profile(
        organism="hsapiens",  # [cite: 231]
        query=query_genes,
        background=bg_genes,  # [cite: 59]
        no_evidences=False,
    )

    if not res.empty:
        res["module_id"] = mod_id
        all_enrichments.append(res)

if all_enrichments:
    pd.concat(all_enrichments).to_csv(snakemake.output["out"], sep="\t", index=False)
else:
    pd.DataFrame().to_csv(snakemake.output["out"], sep="\t", index=False)
