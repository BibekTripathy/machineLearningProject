from pathlib import Path
import pandas as pd

# Load final results
bg_mod = pd.read_csv(snakemake.input["bg_mod"], sep="\t")
dis_net = pd.read_csv(snakemake.input["dis_net"], sep="\t")
enrich = pd.read_csv(snakemake.input["enrich"], sep="\t")

# Basic stats
n_modules = bg_mod["module_id"].nunique()
n_genes = bg_mod["gene_symbol"].nunique()
n_edges = len(dis_net)

# Style
style = """
<style>
    body { font-family: sans-serif; margin: 2rem; color: #333; line-height: 1.6; }
    h1 { color: #2c3e50; border-bottom: 2px solid #eee; padding-bottom: 0.5rem; }
    h2 { color: #34495e; margin-top: 2rem; border-bottom: 1px solid #eee; padding-bottom: 0.2rem;}
    table { border-collapse: collapse; width: 100%; margin-bottom: 1rem; font-size: 0.9rem; }
    th, td { text-align: left; padding: 12px 8px; border-bottom: 1px solid #ddd; }
    th { background-color: #f8f9fa; font-weight: 600; color: #495057; }
    tr:hover { background-color: #f1f1f1; }
    .metric { font-weight: bold; color: #27ae60; }
    .container { max_width: 1200px; margin: 0 auto; }
</style>
"""

# Module Summary Table
# Group by module to get list of genes and size
module_sizes = bg_mod.groupby("module_id")["gene_symbol"].apply(list).reset_index()
module_sizes["Size"] = module_sizes["gene_symbol"].apply(len)
# Create a string of top 10 genes for display
module_sizes["Top Genes"] = module_sizes["gene_symbol"].apply(lambda x: ", ".join(x[:10]) + ("..." if len(x)>10 else ""))
# Rename for display
module_sizes = module_sizes.rename(columns={"module_id": "Module ID"})
module_summary_html = module_sizes[["Module ID", "Size", "Top Genes"]].to_html(index=False, classes="table", border=0)

# Enrichment Table (Top 5 per module)
# We sort by p-value first to ensure we get the most significant ones
enrich_sorted = enrich.sort_values(["module_id", "p_value"])
top_enrich = enrich_sorted.groupby("module_id").head(5)
# Select and rename columns for clarity
display_cols = ["module_id", "source", "native", "name", "p_value"]
rename_map = {
    "module_id": "Module ID",
    "source": "Source",
    "native": "Term ID",
    "name": "Term Name",
    "p_value": "P-Value"
}
enrich_html = top_enrich[display_cols].rename(columns=rename_map).to_html(index=False, classes="table", border=0, float_format="%.2e")

html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Phase 4 Final Report</title>
    {style}
</head>
<body>
    <div class="container">
        <h1>Phase 4 Final Report: Cancer Stratification Analysis</h1>
        
        <h2>1. Executive Summary</h2>
        <p>
            This report summarizes the results of the network propagation and community detection pipeline.
        </p>
        <ul>
            <li>Total Genes Analyzed (in modules): <span class="metric">{n_genes}</span></li>
            <li>Identified Modules: <span class="metric">{n_modules}</span></li>
            <li>Disease-Specific Interactions (Edges): <span class="metric">{n_edges}</span></li>
        </ul>

        <h2>2. Module Composition</h2>
        <p>
            The network was partitioned into distinct modules (communities) using the Leiden algorithm.
            Below is the composition of each identified module.
        </p>
        {module_summary_html}

        <h2>3. Functional Enrichment (Top Pathways)</h2>
        <p>
            Functional enrichment analysis (using g:Profiler) reveals the biological processes associated with each module.
            Showing the top 5 most significant terms (p-value < 0.05) per module.
        </p>
        {enrich_html}

        <h2>4. Next Steps (Phase 6)</h2>
        <p>
            These results constitute the final deliverable for Phase 4.
            The identified modules and their enriched pathways provide a basis for stratified cancer analysis.
        </p>
    </div>
</body>
</html>
"""
Path(snakemake.output["html"]).write_text(html_content)