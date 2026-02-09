# WGCNA Module

## Purpose

Weighted Gene Co-expression Network Analysis (WGCNA) to identify modules of co-expressed genes and correlate them with sample traits including ancestry proportions, clinical variables, and molecular features.

## Processes

| Process | Description |
|---------|-------------|
| `WGCNA_ANALYSIS` | Construct signed co-expression network, detect gene modules via dynamic tree cutting, compute module eigengenes, test module-trait correlations (ancestry, relapse, subtype, etc.), and identify hub genes |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| WGCNA | (Bioc 3.19) | [horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) |

## Inputs

- `path normalized_counts` -- Normalized expression matrix (VST or similar)
- `path metadata` -- Sample metadata with clinical and molecular traits
- `path ancestry_proportions` -- Continuous ancestry proportions for trait correlation

## Outputs

- `wgcna_results/module_eigengenes.tsv` -- Module eigengene values per sample
- `wgcna_results/module_membership.tsv` -- Gene-to-module assignment with membership scores
- `wgcna_results/module_trait_cor.tsv` -- Module-trait correlation matrix with p-values
- `wgcna_results/hub_genes.tsv` -- Top hub genes per module
- `wgcna_plots/` -- Scale-free topology fit, dendrogram, module-trait heatmap

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `wgcna_min_module_size` | `30` | Minimum number of genes per module |
| `wgcna_merge_cut_height` | `0.25` | Height cut for merging similar modules |
| `wgcna_soft_power` | `null` | Soft-thresholding power (auto-detected if null) |
| `wgcna_network_type` | `signed` | Network type: signed or unsigned |
| `wgcna_tom_type` | `signed` | Topological Overlap Matrix type |
