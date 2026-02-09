#!/usr/bin/env python3
"""
Generate integrated HTML report summarizing all pipeline outputs.
Compiles key findings from each module into a single navigable report.
"""

import argparse
import os
import glob
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(description='Generate integrated pipeline report')
    parser.add_argument('--results-dir', default='.')
    parser.add_argument('--output-html', default='integrated_report.html')
    parser.add_argument('--output-tables', default='summary_tables')
    return parser.parse_args()

def find_summary_files(results_dir):
    """Find summary/key output files from each module."""
    summaries = {}
    patterns = {
        'de_summary': '**/de_summary.tsv',
        'limma_summary': '**/limma_summary.tsv',
        'pathway_gsea': '**/gsea/*.tsv',
        'immune_deconv': '**/deconvolution_all.tsv',
        'ancestry_proportions': '**/ancestry_proportions.tsv',
        'ancestry_categories': '**/ancestry_categories.tsv',
        'fusion_summary': '**/fusion_summary.tsv',
        'hla_genotypes': '**/hla_genotypes.tsv',
        'neoantigen_burden': '**/neoantigen_burden.tsv',
        'tmb_scores': '**/tmb_scores.tsv',
        'drug_scores': '**/drug_sensitivity_scores.tsv',
        'tcr_summary': '**/tcr_repertoire_summary.tsv',
        'wgcna_hub_genes': '**/hub_genes.tsv',
        'sensitivity': '**/timepoint_analysis.tsv',
    }

    for key, pattern in patterns.items():
        matches = glob.glob(os.path.join(results_dir, pattern), recursive=True)
        if matches:
            summaries[key] = matches[0]

    return summaries

def read_tsv_as_html_table(filepath, max_rows=50):
    """Read a TSV file and convert to HTML table."""
    try:
        with open(filepath) as f:
            lines = f.readlines()
        if not lines:
            return '<p>Empty file</p>'

        header = lines[0].strip().split('\t')
        rows = [line.strip().split('\t') for line in lines[1:max_rows+1]]

        html = '<table class="data-table">\n<thead><tr>'
        html += ''.join(f'<th>{h}</th>' for h in header)
        html += '</tr></thead>\n<tbody>\n'
        for row in rows:
            html += '<tr>' + ''.join(f'<td>{c}</td>' for c in row) + '</tr>\n'
        html += '</tbody></table>\n'

        if len(lines) > max_rows + 1:
            html += f'<p class="note">Showing {max_rows} of {len(lines)-1} rows</p>\n'

        return html
    except Exception as e:
        return f'<p class="error">Error reading {filepath}: {e}</p>'

def generate_html(summaries, output_html):
    """Generate the integrated HTML report."""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Cancer RNA-Seq Pipeline Report</title>
<style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
           max-width: 1200px; margin: 0 auto; padding: 20px; background: #f5f5f5; }}
    h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
    h2 {{ color: #34495e; border-bottom: 1px solid #bdc3c7; padding-bottom: 5px; margin-top: 30px; }}
    h3 {{ color: #7f8c8d; }}
    .section {{ background: white; padding: 20px; margin: 15px 0; border-radius: 5px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
    .data-table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
    .data-table th, .data-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
    .data-table th {{ background: #3498db; color: white; }}
    .data-table tr:nth-child(even) {{ background: #f2f2f2; }}
    .data-table tr:hover {{ background: #e6f3ff; }}
    .note {{ color: #7f8c8d; font-style: italic; }}
    .error {{ color: #e74c3c; }}
    .nav {{ position: fixed; top: 0; left: 0; width: 200px; height: 100vh; background: #2c3e50;
            color: white; padding: 20px; overflow-y: auto; }}
    .nav a {{ color: #ecf0f1; text-decoration: none; display: block; padding: 5px 0; font-size: 14px; }}
    .nav a:hover {{ color: #3498db; }}
    .content {{ margin-left: 240px; }}
    .badge {{ display: inline-block; padding: 2px 8px; border-radius: 3px; font-size: 12px;
              font-weight: bold; margin: 2px; }}
    .badge-success {{ background: #27ae60; color: white; }}
    .badge-warning {{ background: #f39c12; color: white; }}
    .badge-info {{ background: #3498db; color: white; }}
</style>
</head>
<body>
<nav class="nav">
    <h3>Pipeline Report</h3>
    <a href="#overview">Overview</a>
    <a href="#ancestry">Ancestry</a>
    <a href="#de">Differential Expression</a>
    <a href="#fusions">Fusions</a>
    <a href="#immune">Immune Analysis</a>
    <a href="#hla">HLA Typing</a>
    <a href="#neoantigen">Neoantigens</a>
    <a href="#tmb">TMB</a>
    <a href="#wgcna">WGCNA</a>
    <a href="#tcr">TCR Repertoire</a>
    <a href="#pharma">Pharmacogenomics</a>
    <a href="#sensitivity">Sensitivity</a>
</nav>

<div class="content">
<h1>Cancer RNA-Seq Analysis Pipeline Report</h1>
<p class="note">Generated: {timestamp}</p>

<div class="section" id="overview">
<h2>Pipeline Overview</h2>
<p>Modules completed: <span class="badge badge-success">{len(summaries)} of 14</span></p>
<ul>
"""

    module_status = {
        'de_summary': 'Differential Expression',
        'pathway_gsea': 'Pathway Enrichment',
        'immune_deconv': 'Immune Deconvolution',
        'ancestry_proportions': 'Ancestry Inference',
        'fusion_summary': 'Fusion Detection',
        'hla_genotypes': 'HLA Typing',
        'neoantigen_burden': 'Neoantigen Prediction',
        'tmb_scores': 'TMB Estimation',
        'drug_scores': 'Pharmacogenomics',
        'tcr_summary': 'TCR Repertoire',
        'wgcna_hub_genes': 'WGCNA',
        'sensitivity': 'Sensitivity Analysis',
    }

    for key, label in module_status.items():
        status = 'badge-success' if key in summaries else 'badge-warning'
        status_text = 'Complete' if key in summaries else 'Not run'
        html += f'<li>{label}: <span class="badge {status}">{status_text}</span></li>\n'

    html += '</ul></div>\n'

    # Add each module section
    sections = [
        ('ancestry', 'Ancestry Inference', ['ancestry_proportions', 'ancestry_categories']),
        ('de', 'Differential Expression', ['de_summary', 'limma_summary']),
        ('fusions', 'Fusion Detection', ['fusion_summary']),
        ('immune', 'Immune Analysis', ['immune_deconv']),
        ('hla', 'HLA Typing', ['hla_genotypes']),
        ('neoantigen', 'Neoantigen Prediction', ['neoantigen_burden']),
        ('tmb', 'Tumor Mutational Burden', ['tmb_scores']),
        ('wgcna', 'WGCNA Co-expression', ['wgcna_hub_genes']),
        ('tcr', 'TCR/BCR Repertoire', ['tcr_summary']),
        ('pharma', 'Pharmacogenomics', ['drug_scores']),
        ('sensitivity', 'Sensitivity Analysis', ['sensitivity']),
    ]

    for section_id, title, keys in sections:
        html += f'<div class="section" id="{section_id}">\n<h2>{title}</h2>\n'
        found = False
        for key in keys:
            if key in summaries:
                html += f'<h3>{key.replace("_", " ").title()}</h3>\n'
                html += read_tsv_as_html_table(summaries[key])
                found = True
        if not found:
            html += '<p class="note">No results available for this module.</p>\n'
        html += '</div>\n'

    html += '</div></body></html>'

    with open(output_html, 'w') as f:
        f.write(html)

    print(f"Report generated: {output_html}")

def main():
    args = parse_args()
    os.makedirs(args.output_tables, exist_ok=True)
    summaries = find_summary_files(args.results_dir)
    print(f"Found {len(summaries)} summary files")
    generate_html(summaries, args.output_html)

if __name__ == '__main__':
    main()
