# Ancestry Module

## Purpose

Infer genetic ancestry from bulk RNA-seq BAMs using GRAF-anc ancestry-informative SNPs (282K sites) and perform sample QC/relatedness checking with Somalier. Produces both continuous ancestry proportions and discrete ancestry categories for use as covariates in downstream analyses.

## Processes

| Process | Description |
|---------|-------------|
| `EXTRACT_GRAF_SNPS` | Call genotypes at 282K GRAF-anc SNP positions using bcftools mpileup/call; extract allele counts (typical yield: 20K-50K SNPs from RNA-seq) |
| `MERGE_GRAF_VCFS` | Merge per-sample VCFs into a multi-sample VCF with bcftools merge |
| `GRAFANC_RUN` | Run GRAF-anc for continental (8 groups) and subcontinental (38 groups) ancestry assignment |
| `ANCESTRY_INFERENCE` | Supervised classification producing continuous proportions and discrete categories via PCA + reference panel projection |
| `SOMALIER_EXTRACT` | Extract genotypes at ~17K coding-region SNPs per sample for QC |
| `SOMALIER_RELATE` | Pairwise relatedness and sample-swap detection across all samples |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| bcftools | 1.19 | [samtools.github.io/bcftools](https://samtools.github.io/bcftools/) |
| GRAF-anc | -- | [GRAF-pop/GRAF-anc (NCBI)](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/Software.cgi) |
| Somalier | 0.2.19 | [github.com/brentp/somalier](https://github.com/brentp/somalier) |
| scikit-learn | -- | [scikit-learn.org](https://scikit-learn.org/) |

## Inputs

- `tuple val(meta), path(bam), path(bai)` -- Sorted, indexed BAM per sample
- `path graf_snp_positions` -- GRAF-anc 282K SNP positions (BED or AncSnpPopAFs.txt)
- `path fasta` -- Reference genome FASTA
- `path grafanc_data` -- GRAF-anc data directory
- `path reference_panel` / `path reference_labels` -- Reference panel for supervised classification
- `path sites_vcf` -- Somalier sites VCF (~17K SNPs)

## Outputs

- `*.graf_snps.vcf.gz` / `*.graf_allele_counts.tsv` -- Per-sample genotypes and allele counts
- `grafanc_results.txt` / `grafanc_ancestry_summary.tsv` -- GRAF-anc ancestry assignments
- `ancestry_proportions.tsv` -- Continuous ancestry proportions (5 components)
- `ancestry_categories.tsv` -- Discrete GRAF categories (EUR, AFR_AM, LA1, LA2, EAS, SAS)
- `somalier.samples.tsv` / `somalier.pairs.tsv` / `somalier.html` -- Relatedness and QC

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `graf_ancestry_categories` | `EUR,AFR_AM,LA1,LA2,EAS,SAS` | Ancestry categories for classification |
| `ancestry_continuous_vars` | `pct_african,pct_amerindigenous,...` | Continuous proportion variable names |
| `somalier_sites` | `null` | Path to Somalier sites VCF |
