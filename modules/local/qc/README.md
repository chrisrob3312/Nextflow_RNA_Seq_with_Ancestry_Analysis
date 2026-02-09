# QC Module

## Purpose

Comprehensive BAM-level quality control and aggregated reporting for RNA-seq samples. Generates alignment statistics, RNA-seq-specific metrics, and a unified MultiQC report.

## Processes

| Process | Description |
|---------|-------------|
| `SAMTOOLS_FLAGSTAT` | Alignment flag summary (mapped, paired, duplicates) |
| `SAMTOOLS_IDXSTATS` | Per-chromosome read counts from BAM index |
| `SAMTOOLS_STATS` | Detailed alignment statistics (insert size, base quality, coverage) |
| `RSEQC_BAMSTAT` | BAM-level QC metrics via RSeQC `bam_stat.py` |
| `RSEQC_READDISTRIBUTION` | Read distribution across genomic features (CDS, UTR, intron, intergenic) |
| `RSEQC_INFEREXPERIMENT` | Infer library strandedness from splice-aware alignments |
| `PICARD_COLLECTRNASEQMETRICS` | RNA-seq-specific metrics: rRNA rate, coding/UTR/intergenic fractions, 5'-to-3' bias |
| `MULTIQC` | Aggregate all QC outputs into a single interactive HTML report |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| Samtools | 1.19 | [htslib.org](https://www.htslib.org/doc/samtools.html) |
| RSeQC | 5.0.3 | [rseqc.sourceforge.net](https://rseqc.sourceforge.net/) |
| Picard | 3.1.1 | [broadinstitute.github.io/picard](https://broadinstitute.github.io/picard/) |
| MultiQC | 1.21 | [multiqc.info](https://multiqc.info/) |

## Inputs

- `tuple val(meta), path(bam), path(bai)` -- Sorted, indexed BAM with sample metadata
- `path gene_bed` -- BED12 gene model file (RSeQC read distribution and infer experiment)
- `path fasta` / `path gtf` -- Reference genome and annotation (Picard)

## Outputs

- `*.flagstat`, `*.idxstats`, `*.stats` -- Samtools metric files
- `*.bam_stat.txt`, `*.read_distribution.txt`, `*.infer_experiment.txt` -- RSeQC reports
- `*.rna_metrics` -- Picard CollectRnaSeqMetrics output
- `*multiqc_report.html` -- Aggregated MultiQC HTML report

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `skip_rseqc` | `false` | Skip RSeQC processes |
| `skip_picard` | `false` | Skip Picard CollectRnaSeqMetrics |
| `multiqc_config` | `null` | Custom MultiQC configuration YAML |
