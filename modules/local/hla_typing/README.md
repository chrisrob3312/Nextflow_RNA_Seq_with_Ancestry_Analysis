# HLA Typing Module

## Purpose

HLA class I and class II allele typing from RNA-seq BAMs using arcasHLA and OptiType. HLA types are required for downstream neoantigen prediction.

## Processes

| Process | Description |
|---------|-------------|
| `ARCASHLA_EXTRACT` | Extract HLA-mapped reads from BAM into paired FASTQ files |
| `ARCASHLA_GENOTYPE` | Genotype HLA alleles (A, B, C, DPB1, DQB1, DQA1, DRB1) to 4-digit resolution |
| `ARCASHLA_MERGE` | Merge per-sample HLA genotype JSONs into a single cohort TSV |
| `OPTITYPE` | Class I HLA typing (A, B, C) using integer linear programming on HLA-region reads |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| arcasHLA | 0.6.0 | [github.com/RabadanLab/arcasHLA](https://github.com/RabadanLab/arcasHLA) |
| OptiType | 1.3.5 | [github.com/FRED-2/OptiType](https://github.com/FRED-2/OptiType) |

## Inputs

- `tuple val(meta), path(bam), path(bai)` -- Sorted, indexed BAM per sample

## Outputs

- `*.extracted.1.fq.gz` / `*.extracted.2.fq.gz` -- HLA-region reads (arcasHLA)
- `*.genotype.json` -- Per-sample HLA genotype (arcasHLA)
- `hla_genotypes.tsv` -- Merged cohort-level HLA genotype table
- `*_optitype/` -- OptiType results including `*_result.tsv` with class I calls

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `hla_tool` | `arcashla` | Tool: `arcashla`, `optitype`, or `both` |
| `arcashla_genes` | `A,B,C,DPB1,DQB1,DQA1,DRB1` | HLA genes to genotype |

## Notes

- arcasHLA types both class I and class II alleles from RNA-seq
- OptiType is class I only but often considered a gold standard for class I
- Running both tools allows cross-validation of class I calls
