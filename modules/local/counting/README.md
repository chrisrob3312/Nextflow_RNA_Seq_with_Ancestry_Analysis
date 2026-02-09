# Counting Module

## Purpose

Gene-level read quantification from BAM files using Subread featureCounts, followed by count matrix merging, normalization (VST/rlog/TMM/TPM), and optional ComBat-seq batch correction.

## Processes

| Process | Description |
|---------|-------------|
| `SUBREAD_FEATURECOUNTS` | Per-sample gene-level read counting against GTF annotation; supports paired-end, strand-specific counting |
| `MERGE_COUNTS` | Combine per-sample count files into a single raw count matrix; filter low-expression genes |
| `NORMALIZE_COUNTS` | Apply variance-stabilizing transform (VST), rlog, TMM, or TPM normalization; generate QC plots |
| `BATCH_CORRECTION` | ComBat-seq batch correction on raw counts preserving count structure for downstream DE |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| Subread (featureCounts) | 2.0.6 | [subread.sourceforge.net](https://subread.sourceforge.net/) |
| DESeq2 | (Bioc 3.19) | [bioconductor.org/packages/DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) |
| sva (ComBat-seq) | (Bioc 3.19) | [bioconductor.org/packages/sva](https://bioconductor.org/packages/release/bioc/html/sva.html) |

## Inputs

- `tuple val(meta), path(bam), path(bai)` -- Sorted, indexed BAM per sample
- `path gtf` -- Gene annotation GTF
- `path metadata` -- Sample metadata CSV/TSV (for normalization and batch correction)

## Outputs

- `*.featureCounts.txt` / `*.featureCounts.txt.summary` -- Per-sample counts and assignment stats
- `raw_count_matrix.tsv` -- Merged gene-by-sample raw count matrix
- `normalized_counts_*.tsv` / `normalized_counts_tpm.tsv` -- Normalized expression matrices
- `batch_corrected_counts.tsv` -- Batch-corrected count matrix
- `normalization_qc/` / `batch_correction_qc/` -- Diagnostic plots (PCA, density, boxplots)

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fc_strandedness` | `2` | Strand-specificity: 0=unstranded, 1=stranded, 2=reverse |
| `fc_count_type` | `exon` | GTF feature type to count |
| `fc_group_features` | `gene_id` | GTF attribute for grouping features |
| `min_gene_counts` | `10` | Minimum total counts to retain a gene |
| `min_samples_expressing` | `3` | Minimum samples with non-zero counts |
| `normalization_method` | `vst` | Normalization: vst, rlog, tmm, or tpm |
| `batch_correction` | `combat_seq` | Batch correction method: combat_seq, sva, or none |
| `batch_variable` | `batch` | Metadata column identifying batch |
