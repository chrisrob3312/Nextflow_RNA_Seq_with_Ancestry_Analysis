# Pathway Enrichment Module

## Purpose

Functional enrichment analysis of differentially expressed genes using over-representation analysis (ORA), gene set enrichment analysis (GSEA), and per-sample pathway activity scoring (GSVA). Queries GO, KEGG, Reactome, Hallmark, and ImmuneSigDB collections from MSigDB.

## Processes

| Process | Description |
|---------|-------------|
| `PATHWAY_ENRICHMENT` | Run clusterProfiler ORA on significant DE gene lists and fgsea-based GSEA on ranked gene lists across all configured databases; generates dotplots, enrichment maps, and ridgeplots |
| `GSVA_ANALYSIS` | Compute per-sample gene set variation analysis (GSVA) scores for all configured MSigDB collections; test for differential pathway activity across groups |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| clusterProfiler | (Bioc 3.19) | [bioconductor.org/packages/clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) |
| fgsea | (Bioc 3.19) | [bioconductor.org/packages/fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) |
| GSVA | (Bioc 3.19) | [bioconductor.org/packages/GSVA](https://bioconductor.org/packages/release/bioc/html/GSVA.html) |
| msigdbr | (CRAN) | [cran.r-project.org/package=msigdbr](https://cran.r-project.org/package=msigdbr) |

## Inputs

- `path de_results` -- Directory of DE result tables (one per contrast)
- `path count_matrix` -- Raw count matrix (for ORA background)
- `path normalized_counts` -- Normalized expression matrix (for GSVA)
- `path metadata` -- Sample metadata for GSVA group comparisons

## Outputs

- `pathway_results/ora/` -- ORA results per contrast and database
- `pathway_results/gsea/` -- GSEA results with NES, p-values, leading edge genes
- `pathway_plots/` -- Dotplots, enrichment maps, ridgeplots
- `gsva_scores.tsv` -- Sample-by-pathway activity score matrix
- `gsva_results/` / `gsva_plots/` -- Differential activity tests and plots

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pathway_databases` | `GO_BP,GO_MF,KEGG,REACTOME,HALLMARK,IMMUNESIGDB` | MSigDB collections to query |
| `msigdb_species` | `Homo sapiens` | Species for msigdbr |
| `gsea_min_size` | `15` | Minimum gene set size |
| `gsea_max_size` | `500` | Maximum gene set size |
| `padj_threshold` | `0.05` | Significance cutoff for ORA input |
| `run_gsva` | `true` | Enable GSVA analysis |
