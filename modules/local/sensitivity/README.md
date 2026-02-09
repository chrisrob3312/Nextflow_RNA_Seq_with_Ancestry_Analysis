# Sensitivity Analysis Module

## Purpose

Assess the robustness of differential expression results to analytical choices including timepoint stratification, ancestry proportion thresholds, covariate inclusion/exclusion, and model specification variations.

## Processes

| Process | Description |
|---------|-------------|
| `SENSITIVITY_ANALYSIS` | Run systematic sensitivity assessments: (1) timepoint-stratified DE (diagnostic vs. relapse), (2) ancestry proportion threshold sweeps, (3) covariate leave-one-out impact analysis, (4) model comparison (DESeq2 vs. limma concordance under different specifications); generates summary tables and visualizations |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| R | (Bioc 3.19) | Custom analysis script using DESeq2/limma |

## Inputs

- `path count_matrix` -- Raw count matrix
- `path metadata` -- Sample metadata with all covariates
- `path ancestry_proportions` -- Continuous ancestry proportions
- `path de_results` -- Primary DE results for comparison baseline

## Outputs

- `sensitivity_results/timepoint_analysis.tsv` -- DE results stratified by timepoint
- `sensitivity_results/ancestry_sensitivity.tsv` -- Effect of ancestry threshold variations
- `sensitivity_results/covariate_impact.tsv` -- Covariate leave-one-out analysis
- `sensitivity_results/model_comparison.tsv` -- Cross-model concordance metrics
- `sensitivity_plots/` -- Upset plots, Jaccard similarity heatmaps, forest plots

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sensitivity_variables` | `timepoint,ancestry_proportion` | Variables to assess |
| `sensitivity_timepoints` | `diagnostic,relapse` | Timepoints for stratified analysis |
| `bootstrap_iterations` | `1000` | Bootstrap iterations for stability assessment |
| `de_covariates` | `batch,sex,age,...` | Baseline covariates for leave-one-out |
