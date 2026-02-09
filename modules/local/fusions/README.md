# Fusions Module

## Purpose

Gene fusion detection from RNA-seq using two independent callers (FusionCatcher and Arriba) with multi-caller merging to identify high-confidence fusions. Includes visualization of fusion breakpoints.

## Processes

| Process | Description |
|---------|-------------|
| `FUSIONCATCHER` | Detect fusions via split-read and paired-end mapping; converts BAM to FASTQ internally; uses built-in database of known oncogenic fusions |
| `ARRIBA` | Detect fusions from STAR-aligned BAMs using split-read and discordant-pair evidence; supports blacklist and known-fusion filtering |
| `ARRIBA_VISUALIZATION` | Generate publication-quality PDF plots of fusion breakpoints with protein domain annotations |
| `MERGE_FUSIONS` | Merge and prioritize calls from both callers; flag fusions detected by multiple tools as high-confidence |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| FusionCatcher | 1.33 | [github.com/ndaniel/fusioncatcher](https://github.com/ndaniel/fusioncatcher) |
| Arriba | 2.4.0 | [github.com/suhrig/arriba](https://github.com/suhrig/arriba) |
| Samtools | 1.19 | [htslib.org](https://www.htslib.org/) (BAM-to-FASTQ for FusionCatcher) |

## Inputs

- `tuple val(meta), path(bam), path(bai)` -- STAR-aligned sorted BAM per sample
- `path fusioncatcher_data` -- FusionCatcher reference database directory
- `path fasta` / `path gtf` -- Reference genome and annotation (Arriba)
- `path blacklist` -- Arriba blacklist of recurrent artifacts
- `path known_fusions` -- Known fusion database (Arriba)
- `path protein_domains` -- Protein domain annotation (Arriba visualization)

## Outputs

- `*_fusioncatcher/` -- FusionCatcher results including `final-list_candidate-fusion-genes.txt`
- `*.arriba.fusions.tsv` -- Arriba fusion calls
- `*.arriba.pdf` -- Arriba fusion visualization PDF
- `merged_fusions.tsv` -- Combined calls from both tools
- `high_confidence_fusions.tsv` -- Fusions called by both tools
- `fusion_summary.tsv` -- Per-sample fusion counts and recurrence

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fusion_tool` | `both` | Run `fusioncatcher`, `arriba`, or `both` |
| `min_fusion_reads` | `3` | Minimum supporting reads to report a fusion |
| `fusioncatcher_data` | `null` | Path to FusionCatcher database |
| `arriba_blacklist` | `null` | Arriba blacklist file |
| `arriba_known_fusions` | `null` | Known fusions for Arriba |
