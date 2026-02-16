from pathlib import Path
import pandas as pd
import json
from datetime import datetime

# Load final results
bg_mod = pd.read_csv(snakemake.input["bg_mod"], sep="\t")
dis_net = pd.read_csv(snakemake.input["dis_net"], sep="\t")
enrich = pd.read_csv(snakemake.input["enrich"], sep="\t")
comparison = pd.read_csv(snakemake.input["comparison"], sep="\t")

# Basic stats
n_modules = bg_mod["module_id"].nunique()
n_genes = bg_mod["gene_symbol"].nunique()
n_edges = len(dis_net)

# Calculate module quality metrics for ML
module_stats = bg_mod.groupby("module_id").agg({
    "gene_symbol": "count",
    "module_density": "mean"
}).reset_index()
module_stats.columns = ["Module ID", "Size", "Avg Density"]

# Comparison Metrics Table
comparison_display = comparison[["Metric", "Value"]].copy()
comparison_html = comparison_display.to_html(index=False, classes="table", border=0)

# Style
style = """
<style>
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; margin: 2rem; color: #333; line-height: 1.6; }
    h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 0.5rem; }
    h2 { color: #34495e; margin-top: 2rem; border-bottom: 1px solid #eee; padding-bottom: 0.3rem;}
    h3 { color: #7f8c8d; margin-top: 1.5rem; }
    table { border-collapse: collapse; width: 100%; margin-bottom: 1rem; font-size: 0.85rem; }
    th, td { text-align: left; padding: 10px 8px; border-bottom: 1px solid #ddd; }
    th { background-color: #3498db; color: white; font-weight: 600; }
    tr:hover { background-color: #f8f9fa; }
    .metric { font-weight: bold; color: #27ae60; font-size: 1.1em; }
    .container { max-width: 1400px; margin: 0 auto; }
    .info-box { background: #f8f9fa; border-left: 4px solid #3498db; padding: 1rem; margin: 1rem 0; }
    .ml-features { background: #fff3cd; border-left: 4px solid #ffc107; padding: 1rem; margin: 1rem 0; }
    code { background: #f4f4f4; padding: 2px 6px; border-radius: 3px; font-family: monospace; }
</style>
"""

# Module Summary Table with quality metrics
module_summary = bg_mod.groupby("module_id").agg({
    "gene_symbol": ["count", lambda x: list(x)],
    "module_density": "mean"
}).reset_index()
module_summary.columns = ["Module ID", "Size", "Genes", "Avg Density"]
module_summary["Top Genes"] = module_summary["Genes"].apply(lambda x: ", ".join(x[:8]) + ("..." if len(x)>8 else ""))
module_summary["Density"] = module_summary["Avg Density"].apply(lambda x: f"{x:.4f}")

# Create ML feature summary
ml_features = module_summary[["Module ID", "Size", "Density"]].copy()
ml_features["Feature Type"] = "module_membership"
ml_features["Description"] = ml_features.apply(lambda x: f"Binary feature: gene in {x['Module ID']} (size={x['Size']}, density={x['Density']})", axis=1)
ml_features_html = ml_features[["Module ID", "Size", "Density", "Description"]].to_html(index=False, classes="table", border=0)

module_summary_html = module_summary[["Module ID", "Size", "Density", "Top Genes"]].to_html(index=False, classes="table", border=0)

# Enrichment Table (Top 5 per module)
if not enrich.empty and "module_id" in enrich.columns:
    enrich_sorted = enrich.sort_values(["module_id", "p_value"])
    top_enrich = enrich_sorted.groupby("module_id").head(5)
    display_cols = [c for c in ["module_id", "source", "name", "p_value"] if c in top_enrich.columns]
    rename_map = {
        "module_id": "Module ID",
        "source": "Source",
        "name": "Term Name",
        "p_value": "P-Value"
    }
    enrich_html = top_enrich[display_cols].rename(columns=rename_map).to_html(index=False, classes="table", border=0, float_format="%.2e")
else:
    enrich_html = "<p><em>No enrichment results available</em></p>"

# Network statistics
if "weight" in dis_net.columns:
    avg_score = dis_net["weight"].mean()
    max_score = dis_net["weight"].max()
    min_score = dis_net["weight"].min()
else:
    avg_score = max_score = min_score = 0

html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Cancer Network Analysis - ML Training Data Report</title>
    <meta charset="utf-8">
</head>
<body>
    <div class="container">
        <h1>🔬 Cancer Network Analysis Report</h1>
        <p class="info-box">
            <strong>Purpose:</strong> ML training data generation for cancer gene stratification<br>
            <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br>
            <strong>Framework:</strong> Network propagation + community detection (Leiden)
        </p>

        <h2>📊 Executive Summary</h2>
        <table>
            <tr><th>Metric</th><th>Value</th></tr>
            <tr><td>Total Genes in Modules</td><td class="metric">{n_genes}</td></tr>
            <tr><td>Identified Modules</td><td class="metric">{n_modules}</td></tr>
            <tr><td>Disease PPI Edges</td><td class="metric">{n_edges}</td></tr>
            <tr><td>Avg Edge Score</td><td class="metric">{avg_score:.3f}</td></tr>
            <tr><td>Edge Score Range</td><td>{min_score:.3f} - {max_score:.3f}</td></tr>
        </table>

        <h2>📈 ML Feature Summary</h2>
        <div class="ml-features">
            <strong>Generated Features for ML Model:</strong>
            <ul>
                <li><code>module_membership</code>: Binary indicators for each module (n={n_modules} modules)</li>
                <li><code>propagation_score</code>: RWR score for each gene (continuous, 0-0.05)</li>
                <li><code>module_size</code>: Size of module each gene belongs to (integer)</li>
                <li><code>module_density</code>: Connectivity density of module (float, 0-1)</li>
            </ul>
        </div>
        {ml_features_html}

        <h2>🧬 Module Composition (Leiden Clustering)</h2>
        <p>
            Modules identified using the Leiden algorithm on the STRING background network (score ≥ 0.7).
            Each module represents a functional gene set for ML feature extraction.
        </p>
        {module_summary_html}

        <h2>🔬 Quantitative Comparison (Phase 4A)</h2>
        <p>
            Comparison between background network modules and disease-specific PPI modules.
            High ARI/NMI indicates stable module structure across network types.
        </p>
        {comparison_html}

        <h2>🧪 Functional Enrichment (Top Pathways per Module)</h2>
        <p>
            g:Profiler enrichment analysis showing biological processes associated with each module.
            These can be used for feature interpretation and validation.
        </p>
        {enrich_html}

        <h2>📁 Output Files for ML Pipeline</h2>
        <table>
            <tr><th>File</th><th>Purpose</th><th>ML Usage</th></tr>
            <tr>
                <td><code>results/genes/expanded_genes.tsv</code></td>
                <td>RWR propagation scores</td>
                <td>Continuous features (gene importance)</td>
            </tr>
            <tr>
                <td><code>results/modules/final_modules.tsv</code></td>
                <td>Module assignments</td>
                <td>Binary features (module membership)</td>
            </tr>
            <tr>
                <td><code>results/networks/string_bg.tsv</code></td>
                <td>Background PPI network</td>
                <td>Graph features, edge weights</td>
            </tr>
            <tr>
                <td><code>results/enrichment/enrich_cancer_leiden.tsv</code></td>
                <td>Pathway enrichment</td>
                <td>Feature interpretation, validation</td>
            </tr>
        </table>

        <h2>✅ Pipeline Status</h2>
        <div class="info-box">
            <strong>Phase 1:</strong> ✅ Input genes loaded (40 cancer genes)<br>
            <strong>Phase 2A:</strong> ✅ STRING network retrieved (237 edges, score ≥ 0.7)<br>
            <strong>Phase 2B:</strong> ✅ RWR propagation completed<br>
            <strong>Phase 2D:</strong> ✅ Leiden clustering completed ({n_modules} modules)<br>
            <strong>Phase 3:</strong> ✅ Disease-specific PPI retrieved<br>
            <strong>Phase 4A:</strong> ✅ Module comparison completed<br>
            <strong>Phase 4B:</strong> ✅ Functional enrichment completed<br>
            <strong>Phase 4:</strong> ✅ Report generated
        </div>
    </div>
</body>
</html>
"""
Path(snakemake.output["html"]).write_text(html_content)
