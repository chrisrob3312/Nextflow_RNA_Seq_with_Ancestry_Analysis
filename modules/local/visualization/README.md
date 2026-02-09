# Visualization Module

## Purpose

Generate integrated publication-quality figures across all analysis modules and maintain a structured pipeline analysis log documenting each step with timestamps, inputs, and key findings.

## Processes

| Process | Description |
|---------|-------------|
| `VISUALIZATION_SUMMARY` | Create summary figures spanning DE, pathway, immune, ancestry, WGCNA, and pharmacogenomics results; generates both PNG and PDF formats; produces an HTML figure index for browsing |
| `PIPELINE_LOGGER` | Record a structured log entry for each pipeline step with timestamp, description, input files, and key findings |
| `COMPILE_PIPELINE_LOG` | Aggregate all log entries into a comprehensive pipeline analysis log in markdown format with parameter documentation |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| R (ggplot2, etc.) | (Bioc 3.19) | Visualization via `create_visualizations.R` |
| Python | -- | HTML figure index generation via `generate_figure_index.py` |

## Inputs

- `path de_results` -- Differential expression results
- `path pathway_results` -- Pathway enrichment results
- `path immune_results` -- Immune deconvolution results
- `path ancestry_results` -- Ancestry inference results
- `path wgcna_results` -- WGCNA network results
- `path pharma_results` -- Pharmacogenomics results
- `path metadata` -- Sample metadata
- `path log_entries` / `path pipeline_params` -- Log entries and parameter records

## Outputs

- `figures/` -- Organized figure directories: `de/`, `pathway/`, `immune/`, `ancestry/`, `wgcna/`, `pharma/`, `overview/`
- `figure_index.html` -- Browsable HTML index of all generated figures
- `pipeline_analysis_log.md` -- Comprehensive analysis log with parameters, steps, and findings

## Key Parameters

This module uses no module-specific parameters. Figure format defaults to PNG and PDF.
