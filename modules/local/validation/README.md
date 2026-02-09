# Validation Module

## Purpose

Quality control checks and artifact detection across pipeline outputs. Validates DE concordance between methods, verifies sample sex from expression, detects expression outliers, and assesses genomic inflation in DE p-values.

## Processes

| Process | Description |
|---------|-------------|
| `VALIDATE_DE_CONCORDANCE` | Compare DESeq2 and limma-voom results per contrast: rank correlation, overlap of significant genes, log2FC agreement, and discordant gene flagging |
| `VALIDATE_SEX_CHECK` | Verify reported sample sex against XIST expression and Y-chromosome gene expression; flag mismatches indicating potential sample swaps |
| `VALIDATE_EXPRESSION_OUTLIERS` | Detect outlier samples using PCA distance, pairwise correlation, and standard deviation thresholds; flag samples below minimum correlation |
| `VALIDATE_GENOMIC_INFLATION` | Compute genomic inflation factor (lambda) from DE p-value distributions; generate QQ plots to detect systematic bias |
| `COMPILE_VALIDATION_REPORT` | Aggregate all validation findings into a unified markdown report with summary table |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| R | (Bioc 3.19) | Custom validation scripts |

## Inputs

- `path deseq2_results` / `path limma_results` -- DE results from both methods
- `path count_matrix` / `path normalized_counts` -- Raw and normalized counts
- `path metadata` -- Sample metadata with reported sex
- `path de_results` -- DE results with p-values (for inflation check)
- `val contrast_name` -- Contrast identifier

## Outputs

- `concordance_*/concordance_summary.tsv` -- DESeq2 vs. limma agreement metrics
- `sex_check/` -- Sex check results and mismatch flags
- `expression_outliers/` -- Outlier sample list with PCA and correlation metrics
- `inflation_*/` -- Lambda values and QQ plots per contrast
- `validation_report.md` -- Compiled validation report
- `validation_summary.tsv` -- Summary table of all checks

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `outlier_sd_threshold` | `3` | SD threshold for PCA-based outlier detection |
| `min_sample_correlation` | `0.8` | Minimum pairwise correlation to pass QC |
