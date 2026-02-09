# Pharmacogenomics Module

## Purpose

Predict drug sensitivity from RNA-seq expression profiles using oncoPredict with GDSC/CCLE training data, identify druggable targets via DGIdb, and optionally query Connectivity Map (CMap) for drug repurposing candidates. Compares predicted drug response across ancestry and clinical groups.

## Processes

| Process | Description |
|---------|-------------|
| `PHARMACOGENOMICS` | (1) Predict per-sample IC50/AUC drug sensitivity scores using oncoPredict trained on GDSC/CCLE cell line pharmacogenomics, (2) query DGIdb for druggable gene-drug interactions from DE results, (3) optionally query CMap signatures, (4) compare drug sensitivity predictions between ancestry and clinical groups |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| oncoPredict | (Bioc 3.19) | [github.com/danculib/oncoPredict](https://github.com/danculib/oncoPredict) |
| DGIdb | -- | [dgidb.org](https://www.dgidb.org/) |

## Inputs

- `path normalized_counts` -- Normalized expression matrix
- `path de_results` -- DE result tables (for druggable target identification)
- `path metadata` -- Sample metadata
- `path ancestry_proportions` -- Ancestry proportions for group comparisons

## Outputs

- `pharma_results/drug_sensitivity_scores.tsv` -- Per-sample predicted IC50/AUC for drugs
- `pharma_results/druggable_targets.tsv` -- DE genes with known drug interactions
- `pharma_results/dgidb_interactions.tsv` -- Full DGIdb interaction table
- `pharma_results/cmap_connections.tsv` -- CMap connectivity scores (optional)
- `pharma_results/group_comparisons/` -- Drug sensitivity by ancestry/clinical groups
- `pharma_plots/` -- Drug sensitivity heatmaps, waterfall plots, group comparisons

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `drug_response_db` | `GDSC` | Training database: `GDSC`, `CCLE`, `PRISM`, or `all` |
| `dgidb_interactions` | `true` | Query DGIdb for druggable targets |
| `cmap_signatures` | `null` | Path to CMap signature database |
