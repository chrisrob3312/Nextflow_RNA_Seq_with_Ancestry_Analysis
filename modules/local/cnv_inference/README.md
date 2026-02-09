# CNV Inference Module

## Purpose

Infer large-scale copy number variations (CNVs) from RNA-seq expression data using InferCNV. Compares tumor expression profiles against a normal reference group to detect chromosomal gains, losses, and focal events.

## Processes

| Process | Description |
|---------|-------------|
| `INFERCNV` | Run InferCNV to detect CNVs by comparing gene expression across chromosomal positions between tumor and reference samples; applies noise filtering, HMM-based state prediction, and generates heatmap visualizations |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| inferCNV | (Bioc 3.19) | [github.com/broadinstitute/inferCNV](https://github.com/broadinstitute/inferCNV) |

## Inputs

- `path count_matrix` -- Raw gene-by-sample count matrix
- `path gene_order_file` -- Gene genomic position file (gene, chr, start, end)
- `path annotations_file` -- Sample annotations defining tumor vs. reference groups

## Outputs

- `infercnv_output/` -- Full InferCNV results directory
- `infercnv.png` -- Chromosomal expression heatmap with inferred CNV states
- `infercnv.observations.txt` -- Per-gene, per-sample CNV scores

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cnv_tool` | `infercnv` | CNV inference tool selection |
| `cnv_reference_group` | `null` | Normal reference group label in annotations file |
| `cnv_gene_order_file` | `null` | Path to gene position ordering file |

## Notes

- InferCNV is computationally intensive (labeled `process_very_high`)
- Requires a reference group of normal/non-tumor samples for baseline
- Particularly useful for detecting B-ALL-associated aneuploidies (hyperdiploidy, hypodiploidy)
