# Variant Calling Module

## Purpose

SNV/indel calling from RNA-seq BAMs following GATK best practices for RNA-seq, including split-read handling, base quality score recalibration, and hard filtering. Also estimates tumor mutational burden (TMB).

## Processes

| Process | Description |
|---------|-------------|
| `GATK_SPLITNCIGARREADS` | Split reads spanning splice junctions into exon-level segments; reassign mapping qualities |
| `GATK_BASERECALIBRATOR` | Model systematic base quality errors using known SNP sites (dbSNP) |
| `GATK_APPLYBQSR` | Apply recalibrated base quality scores to BAM |
| `GATK_HAPLOTYPECALLER` | Call SNVs and indels via local de novo assembly with soft-clipped bases excluded |
| `GATK_VARIANTFILTRATION` | Hard-filter variants on FS (>30), QD (<2), MQ (<40), and depth |
| `TMB_ESTIMATION` | Calculate tumor mutational burden from filtered variants per sample |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| GATK | 4.5.0.0 | [gatk.broadinstitute.org](https://gatk.broadinstitute.org/hc/en-us) |
| R | (Bioc 3.19) | Used for TMB estimation script |

## Inputs

- `tuple val(meta), path(bam), path(bai)` -- Sorted, indexed BAM per sample
- `path fasta` / `path fasta_fai` -- Reference genome and index
- `path known_snps` / `path known_snps_tbi` -- dbSNP VCF for BQSR
- `path dbsnp` / `path dbsnp_tbi` -- dbSNP for HaplotypeCaller annotation
- `path metadata` -- Sample metadata (for TMB grouping)

## Outputs

- `*.split.bam` -- BAM with split CIGAR reads (intermediate, not published)
- `*.recal.bam` -- Recalibrated BAM (intermediate, not published)
- `*.raw.vcf.gz` -- Raw variant calls
- `*.filtered.vcf.gz` -- Hard-filtered variant calls
- `tmb_scores.tsv` / `tmb_summary.tsv` -- Per-sample TMB estimates

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_base_quality` | `20` | Minimum base quality for variant calling |
| `min_variant_depth` | `10` | Minimum read depth at variant site |
| `min_vaf` | `0.05` | Minimum variant allele fraction for TMB |
| `dbsnp` | `null` | dbSNP VCF path |
