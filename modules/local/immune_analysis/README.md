# Immune Analysis Module

## Purpose

Estimate immune cell composition and tumor purity from RNA-seq expression data using six deconvolution methods and ESTIMATE. Results are corrected for tumor purity when applicable.

## Processes

| Process | Description |
|---------|-------------|
| `IMMUNE_DECONVOLUTION` | Run up to six immune deconvolution methods (CIBERSORTx, xCell, MCP-counter, EPIC, TIMER, quantiseq) via the immunedeconv R framework; optionally correct cell fractions for tumor purity |
| `ESTIMATE_SCORES` | Compute stromal score, immune score, ESTIMATE score, and inferred tumor purity using the ESTIMATE algorithm |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| immunedeconv | (Bioc 3.19) | [github.com/omnideconv/immunedeconv](https://github.com/omnideconv/immunedeconv) |
| ESTIMATE | (Bioc 3.19) | [bioinformatics.mdanderson.org/estimate](https://bioinformatics.mdanderson.org/estimate/) |
| CIBERSORTx | -- | [cibersortx.stanford.edu](https://cibersortx.stanford.edu/) (requires token) |
| xCell | -- | [github.com/dviraran/xCell](https://github.com/dviraran/xCell) |
| MCP-counter | -- | [github.com/ebecht/MCPcounter](https://github.com/ebecht/MCPcounter) |
| EPIC | -- | [github.com/GfellerLab/EPIC](https://github.com/GfellerLab/EPIC) |

## Inputs

- `path normalized_counts` -- Normalized expression matrix (TPM recommended for most methods)
- `path metadata` -- Sample metadata for grouping and visualization

## Outputs

- `immune_results/deconvolution_all.tsv` -- Combined scores from all methods
- `immune_results/cell_fractions.tsv` -- Estimated cell type fractions
- `immune_plots/` -- Box plots, heatmaps, and correlation plots across methods
- `estimate_scores.tsv` -- Stromal, immune, and ESTIMATE scores
- `estimate_purity.tsv` -- Inferred tumor purity per sample
- `estimate_plots/` -- Purity distribution and score visualizations

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `immune_deconv_methods` | `cibersortx,xcell,mcpcounter,estimate,epic,timer` | Methods to run |
| `cibersortx_token` | `null` | CIBERSORTx API token (required for CIBERSORTx) |
| `cibersortx_sigmatrix` | `LM22` | Signature matrix for CIBERSORTx |
| `estimate_platform` | `illumina` | Platform type for ESTIMATE |
| `correct_tumor_purity` | `true` | Adjust immune fractions for tumor purity |
