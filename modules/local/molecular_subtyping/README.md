# Molecular Subtyping Module

## Purpose

Classify samples into cytomolecular subtypes based on gene expression signatures. Default configuration targets pediatric B-ALL with 12+ established subtypes (ETV6-RUNX1, BCR-ABL1, KMT2A-rearranged, hyperdiploid, hypodiploid, etc.) using curated gene signature sets.

## Processes

| Process | Description |
|---------|-------------|
| `MOLECULAR_SUBTYPING` | Classify each sample into a molecular subtype using expression-based signatures; supports consensus classification across multiple gene set scoring methods; outputs subtype assignments and confidence scores |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| R | (Bioc 3.19) | Custom classification script |

## Inputs

- `path normalized_counts` -- Normalized expression matrix
- `path metadata` -- Sample metadata (may include known subtype for validation)

## Outputs

- `subtyping_results/subtypes.tsv` -- Per-sample subtype assignments
- `subtyping_results/classifier_scores.tsv` -- Per-sample signature scores for each subtype
- `subtyping_plots/` -- Heatmaps of signature scores, UMAP/PCA by subtype

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `subtyping_method` | `consensus` | Classification method |

## B-ALL Subtypes

The default gene signatures cover established B-ALL cytomolecular subtypes including:
- ETV6-RUNX1 (t(12;21))
- BCR-ABL1 / Ph-like
- KMT2A-rearranged (MLL)
- High hyperdiploidy (51-67 chromosomes)
- Low hypodiploidy / near-haploid
- TCF3-PBX1 (t(1;19))
- iAMP21
- DUX4-rearranged
- MEF2D-rearranged
- ZNF384-rearranged
- NUTM1-rearranged
- PAX5alt / PAX5 P80R
