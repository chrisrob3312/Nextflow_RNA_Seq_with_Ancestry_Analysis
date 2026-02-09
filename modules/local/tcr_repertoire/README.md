# TCR/BCR Repertoire Module

## Purpose

Reconstruct T-cell receptor (TCR) and B-cell receptor (BCR) repertoires from bulk RNA-seq using TRUST4. Provides clonotype identification, diversity metrics, and cross-sample clonotype tracking.

## Processes

| Process | Description |
|---------|-------------|
| `TRUST4` | Extract and assemble TCR/BCR CDR3 sequences from RNA-seq BAM reads mapping to immune receptor loci; reports clonotype frequency, V/D/J gene usage, and CDR3 amino acid sequences |
| `MERGE_TCR_REPORTS` | Merge per-sample TRUST4 reports into cohort-level summaries; compute diversity indices (Shannon, Simpson, clonality) and track shared clonotypes across samples |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| TRUST4 | 1.1.0 | [github.com/liulab-dfci/TRUST4](https://github.com/liulab-dfci/TRUST4) |

## Inputs

- `tuple val(meta), path(bam), path(bai)` -- Sorted, indexed BAM per sample
- `path fasta` -- Reference genome FASTA (TRUST4 uses human IMGT+C reference internally)

## Outputs

- `*_TRUST4_report.tsv` -- Per-sample clonotype report (CDR3 sequences, frequency, V/D/J usage)
- `*_TRUST4_barcode_report.tsv` -- Barcode-level report
- `*_TRUST4_annot.fa` -- Assembled receptor sequences in FASTA format
- `tcr_repertoire_summary.tsv` -- Cohort-level TCR/BCR summary
- `tcr_diversity_metrics.tsv` -- Shannon entropy, Simpson index, clonality per sample
- `tcr_clonotype_tracking.tsv` -- Shared clonotypes across samples

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tcr_tool` | `trust4` | TCR/BCR repertoire tool |

## Notes

- TRUST4 can reconstruct both TCR (alpha/beta, gamma/delta) and BCR (heavy/light) chains
- Sensitivity depends on immune cell infiltration; bone marrow samples typically have higher yields
- Diversity metrics enable comparison of immune repertoire complexity between groups
