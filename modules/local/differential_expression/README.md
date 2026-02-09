# Differential Expression Module

## Purpose

Differential gene expression analysis using both DESeq2 and limma-voom for concordance. Supports contrasts on continuous ancestry proportions, categorical ancestry (GRAF groups), relapse status, and ADI quartile, including within-cytomolecular-subgroup analyses.

## Processes

| Process | Description |
|---------|-------------|
| `DESEQ2_DE` | Negative binomial GLM-based differential expression via DESeq2; adjusts for covariates (batch, sex, age, blast%, tumor purity, timepoint); runs all contrasts from JSON |
| `LIMMA_VOOM_DE` | Linear modeling on voom-transformed counts via limma; same covariate adjustment and contrasts as DESeq2 for cross-method concordance |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| DESeq2 | (Bioc 3.19) | [bioconductor.org/packages/DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) |
| limma | (Bioc 3.19) | [bioconductor.org/packages/limma](https://bioconductor.org/packages/release/bioc/html/limma.html) |
| edgeR | (Bioc 3.19) | [bioconductor.org/packages/edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) |

## Inputs

- `path count_matrix` -- Raw gene-by-sample count matrix
- `path metadata` -- Sample metadata with covariate columns
- `path ancestry_proportions` -- Continuous ancestry proportions from ancestry module
- `path contrasts_json` -- JSON defining contrasts (variable, reference, comparison, type)

## Outputs

- `deseq2_results/` / `limma_results/` -- Per-contrast result tables (gene, log2FC, padj, etc.)
- `deseq2_plots/` / `limma_plots/` -- Volcano plots, MA plots, heatmaps per contrast
- `deseq2_rds/` / `limma_rds/` -- Serialized R objects for downstream use

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `de_tool` | `both` | Run `deseq2`, `limma`, or `both` |
| `padj_threshold` | `0.05` | Adjusted p-value significance cutoff |
| `lfc_threshold` | `0.585` | log2 fold-change threshold (log2(1.5)) |
| `de_covariates` | `batch,sex,age,blast_percentage,tumor_purity,timepoint` | Covariates in the model formula |
| `de_contrasts` | `null` | Path to JSON file defining contrasts |
| `cytomolecular_subgroups` | `null` | Subgroups for within-group DE analyses |
