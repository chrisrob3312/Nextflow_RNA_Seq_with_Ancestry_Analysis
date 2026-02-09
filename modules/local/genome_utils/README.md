# Genome Utilities Module

## Purpose

Auto-detect genome build (hg19/GRCh37 vs. hg38/GRCh38) from BAM headers and perform coordinate liftover using CrossMap when input BAMs do not match the target build. Ensures all samples are analyzed on a consistent reference genome.

## Processes

| Process | Description |
|---------|-------------|
| `DETECT_GENOME_BUILD` | Parse BAM header to determine genome build from chromosome naming convention (chr prefix) and contig lengths (chr1: 248956422 = hg38, 249250621 = hg19); also detects UCSC vs. Ensembl chromosome style |
| `CROSSMAP_BAM` | Liftover BAM coordinates between genome builds using a chain file (e.g., hg19ToHg38); re-sorts and re-indexes the lifted BAM; reports mapping statistics |
| `VALIDATE_GENOME_BUILDS` | Aggregate build detection results across all samples; warn if mixed builds are detected |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| Samtools | 1.19 | [htslib.org](https://www.htslib.org/) |
| CrossMap | -- | [crossmap.sourceforge.net](http://crossmap.sourceforge.net/) |

## Inputs

- `tuple val(meta), path(bam), path(bai)` -- Sorted, indexed BAM per sample
- `path target_fasta` -- Target genome FASTA (e.g., hg38)
- `path chain_file` -- Liftover chain file (e.g., hg19ToHg38.over.chain.gz)

## Outputs

- `*.genome_build.txt` -- Per-sample detected build, chromosome style, and chr1 length
- `*.liftover.sorted.bam` / `*.liftover.sorted.bam.bai` -- Lifted-over BAM and index
- `*.liftover_stats.txt` -- Liftover statistics (total reads, mapped after liftover, % mapped)
- `genome_build_summary.tsv` -- Cohort-level build summary with mixed-build warnings

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `target_genome_build` | `hg38` | Target build for analysis: `hg38` or `hg19` |
| `auto_detect_build` | `true` | Auto-detect genome build from BAM header |
| `force_liftover` | `false` | Force liftover even if detected build matches target |
| `chain_file` | `null` | Path to liftover chain file |

## Notes

- Detection uses chr1 contig length as the primary discriminator
- Fallback detection checks for GRCh38/hg38/GRCh37/hg19 strings in BAM header
- CrossMap correctly handles spliced RNA-seq alignments during liftover
- Typical liftover success rate is >99% of mapped reads
