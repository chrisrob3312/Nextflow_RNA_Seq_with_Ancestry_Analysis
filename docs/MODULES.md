# Comprehensive Module Reference

Complete technical reference for all processes in the Cancer Bulk RNA-Seq Analysis Pipeline with Ancestry Integration. Organized by subworkflow execution order.

**Pipeline version:** 1.0.0 | **Nextflow:** >=23.04.0 | **DSL:** 2

---

## Table of Contents

1. [Pipeline Architecture](#1-pipeline-architecture)
2. [Subworkflow 1: Preprocessing](#2-subworkflow-1-preprocessing)
   - [Genome Utilities](#21-genome-utilities)
   - [Quality Control](#22-quality-control)
   - [Read Counting and Normalization](#23-read-counting-and-normalization)
3. [Subworkflow 2: Genomics](#3-subworkflow-2-genomics)
   - [Ancestry Inference](#31-ancestry-inference)
   - [Variant Calling and TMB](#32-variant-calling-and-tmb)
   - [HLA Typing](#33-hla-typing)
4. [Subworkflow 3: Expression Analysis](#4-subworkflow-3-expression-analysis)
   - [Differential Expression](#41-differential-expression)
   - [WGCNA Co-expression Networks](#42-wgcna-co-expression-networks)
   - [Pathway Enrichment](#43-pathway-enrichment)
5. [Subworkflow 4: Structural Variants](#5-subworkflow-4-structural-variants)
   - [Differential Splicing](#51-differential-splicing)
   - [Gene Fusion Detection](#52-gene-fusion-detection)
   - [CNV Inference](#53-cnv-inference)
6. [Subworkflow 5: Immunogenomics](#6-subworkflow-5-immunogenomics)
   - [Immune Deconvolution](#61-immune-deconvolution)
   - [Neoantigen Prediction](#62-neoantigen-prediction)
   - [TCR/BCR Repertoire](#63-tcrbcr-repertoire)
7. [Subworkflow 6: Clinical Analysis](#7-subworkflow-6-clinical-analysis)
   - [Sensitivity Analysis](#71-sensitivity-analysis)
   - [Pharmacogenomics](#72-pharmacogenomics)
   - [Molecular Subtyping](#73-molecular-subtyping)
8. [Validation](#8-validation)
9. [Visualization and Reporting](#9-visualization-and-reporting)
10. [Complete Software Reference](#10-complete-software-reference)

---

## 1. Pipeline Architecture

The pipeline accepts STAR-aligned, sorted, duplicate-marked BAM/BAI files and routes them through six subworkflows. Each subworkflow is independently toggleable.

```
BAM/BAI → PREPROCESSING → GENOMICS      → EXPRESSION_ANALYSIS → STRUCTURAL_VARIANTS
                        → IMMUNOGENOMICS → CLINICAL_ANALYSIS   → VALIDATION & REPORTING
```

**Common input format:** All per-sample modules accept `tuple val(meta), path(bam), path(bai)` where `meta` is a Groovy map containing sample_id, batch, sex, age, tumor_purity, cytomolecular_subgroup, relapse_status, adi_quartile, timepoint, disease_stage, and blast_percentage.

---

## 2. Subworkflow 1: Preprocessing

**File:** `subworkflows/local/preprocessing.nf`

Performs genome build detection, optional hg19-to-hg38 liftover, BAM-level QC, gene counting, normalization, and batch correction.

### 2.1 Genome Utilities

**Module:** `modules/local/genome_utils/main.nf`

Auto-detects genome build from BAM headers and performs coordinate liftover when input BAMs do not match the target build, ensuring all samples are analyzed on a consistent reference genome.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `DETECT_GENOME_BUILD` | Parse BAM header to determine genome build from chr1 contig length (248956422 = hg38, 249250621 = hg19) and chromosome naming convention (UCSC `chr` prefix vs Ensembl numeric). Fallback detection checks for assembly identifiers (GRCh38/GRCh37/GCA accessions) in the header. Outputs the detected build as a Nextflow `env(BUILD)` variable for downstream branching. | `process_low` |
| `CROSSMAP_BAM` | Liftover BAM coordinates between builds using a chain file. Handles spliced RNA-seq alignments correctly via CrossMap. Re-sorts, re-indexes, and collects mapping statistics (total reads, mapped after liftover, % mapped). Typical success rate: >99%. | `process_high` |
| `VALIDATE_GENOME_BUILDS` | Aggregate build detection results across all samples. Warns to stderr if mixed builds are detected in the cohort. | `process_low` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `target_genome_build` | `hg38` | Target build for analysis |
| `auto_detect_build` | `true` | Auto-detect genome build from BAM header |
| `force_liftover` | `false` | Force liftover even if detected build matches target |
| `chain_file` | `null` | Path to liftover chain file (e.g., `hg19ToHg38.over.chain.gz`) |

**Outputs:** Per-sample `*.genome_build.txt` (build, chr style, chr1 length), lifted-over BAMs (intermediate), cohort-level `genome_build_summary.tsv`.

**Software:** Samtools 1.19, CrossMap

---

### 2.2 Quality Control

**Module:** `modules/local/qc/main.nf`

Comprehensive BAM-level QC and aggregated reporting. Generates alignment statistics, RNA-seq-specific metrics, and a unified MultiQC report.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `SAMTOOLS_FLAGSTAT` | Alignment flag summary: mapped, paired, properly paired, duplicates, singletons | `process_low` |
| `SAMTOOLS_IDXSTATS` | Per-chromosome read counts from BAM index | `process_low` |
| `SAMTOOLS_STATS` | Comprehensive alignment statistics: insert size distribution, base quality, coverage depth, error rates | `process_low` |
| `RSEQC_BAMSTAT` | BAM-level QC metrics via `bam_stat.py`: uniquely mapped, multi-mapped, non-primary alignments | `process_low` |
| `RSEQC_READDISTRIBUTION` | Read distribution across genomic features: CDS, 5'UTR, 3'UTR, intron, intergenic regions | `process_medium` |
| `RSEQC_INFEREXPERIMENT` | Infer library strandedness from splice-aware alignments; reports fraction of reads in sense/antisense orientation | `process_low` |
| `PICARD_COLLECTRNASEQMETRICS` | RNA-seq-specific metrics: rRNA contamination rate, coding/UTR/intergenic read fractions, 5'-to-3' coverage bias, median CV of gene coverage | `process_medium` |
| `MULTIQC` | Aggregate all QC outputs across all samples into a single interactive HTML report | `process_low` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `skip_rseqc` | `false` | Skip RSeQC processes |
| `skip_picard` | `false` | Skip Picard CollectRnaSeqMetrics |
| `multiqc_config` | `null` | Custom MultiQC configuration YAML |

**Outputs:** Per-sample `*.flagstat`, `*.idxstats`, `*.stats`; per-sample RSeQC `*.bam_stat.txt`, `*.read_distribution.txt`, `*.infer_experiment.txt`; per-sample Picard `*.rna_metrics`; aggregated `multiqc_report.html`.

**Software:** Samtools 1.19, RSeQC 5.0.3, Picard 3.1.1, MultiQC 1.21

---

### 2.3 Read Counting and Normalization

**Module:** `modules/local/counting/main.nf`

Gene-level read quantification, count matrix assembly, normalization, and optional batch correction.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `SUBREAD_FEATURECOUNTS` | Per-sample gene-level read counting against GTF annotation using Subread featureCounts. Supports paired-end reads, strand-specific counting (reversely stranded by default), and configurable feature types. Reports assignment summary statistics. | `process_medium` |
| `MERGE_COUNTS` | Combine per-sample count files into a single gene-by-sample raw count matrix. Filters low-expression genes based on minimum count and minimum sample thresholds. | `process_low` |
| `NORMALIZE_COUNTS` | Apply one of four normalization methods: variance-stabilizing transform (VST), regularized log (rlog), trimmed mean of M-values (TMM), or transcripts per million (TPM). Generates diagnostic QC plots (PCA, density, boxplots). | `process_medium` |
| `BATCH_CORRECTION` | ComBat-seq batch correction on raw counts, preserving the discrete count structure for downstream DE analysis. Alternative: SVA surrogate variable analysis. | `process_medium` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fc_strandedness` | `2` | Strand-specificity: 0=unstranded, 1=stranded, 2=reversely stranded |
| `fc_count_type` | `exon` | GTF feature type to count |
| `fc_group_features` | `gene_id` | GTF attribute for grouping features into genes |
| `fc_extra_attributes` | `gene_name` | Extra GTF attributes to include in output |
| `min_gene_counts` | `10` | Minimum total counts across samples to retain a gene |
| `min_samples_expressing` | `3` | Minimum number of samples with non-zero counts |
| `normalization_method` | `vst` | Normalization: `vst`, `rlog`, `tmm`, or `tpm` |
| `batch_correction` | `combat_seq` | Batch correction: `combat_seq`, `sva`, or `none` |
| `batch_variable` | `batch` | Metadata column identifying batch |

**Outputs:** Per-sample `*.featureCounts.txt` with summary; `raw_count_matrix.tsv`; `normalized_counts_*.tsv`; `batch_corrected_counts.tsv`; diagnostic QC plots.

**Software:** Subread (featureCounts) 2.0.6, DESeq2, edgeR, sva (ComBat-seq) (Bioc 3.19)

---

## 3. Subworkflow 2: Genomics

**File:** `subworkflows/local/genomics.nf`

Ancestry inference from RNA-seq reads, GATK variant calling, HLA typing, and TMB estimation.

### 3.1 Ancestry Inference

**Module:** `modules/local/ancestry/main.nf`

Infer genetic ancestry from bulk RNA-seq BAMs using GRAF-anc ancestry-informative SNPs (282,424 sites) and perform sample QC/relatedness checking with Somalier. Produces continuous ancestry proportions and discrete ancestry categories for use as DE covariates.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `EXTRACT_GRAF_SNPS` | Call genotypes at 282K GRAF-anc SNP positions using `bcftools mpileup` and `bcftools call`. RNA-seq typically covers 20K-50K of these SNPs (those in expressed regions). Outputs per-sample VCF and allele count tables. | `process_medium` |
| `MERGE_GRAF_VCFS` | Merge per-sample VCFs into a multi-sample VCF using `bcftools merge`. Required for cohort-level GRAF-anc analysis. | `process_medium` |
| `GRAFANC_RUN` | Run the GRAF-anc binary on the merged multi-sample VCF. Uses geometric genetic distance from reference populations (not PCA). Assigns continental (8 groups) and subcontinental (38 groups) ancestry via barycentric coordinates. | `process_medium` |
| `ANCESTRY_INFERENCE` | Supervised classification using GRAF-anc output plus allele counts. Produces: (1) continuous proportions for 5 ancestral components (pct_african, pct_amerindigenous, pct_east_asian, pct_south_asian, pct_european), and (2) discrete GRAF categories (EUR, AFR_AM, LA1, LA2, EAS, SAS). Uses PCA + reference panel projection when reference panel is provided. | `process_medium` |
| `SOMALIER_EXTRACT` | Extract genotypes at ~17K coding-region SNPs per sample for downstream QC. Independent of the GRAF-anc workflow. | `process_low` |
| `SOMALIER_RELATE` | Pairwise relatedness estimation and sample-swap detection across all samples. Outputs interactive HTML report, pairs TSV, and samples TSV. | `process_low` |

**Coverage notes:**
- 10K+ GRAF SNPs from RNA-seq: sufficient for continental-level ancestry (EUR/AFR/EAS/SAS/AMR)
- 20K+ GRAF SNPs: recommended for subcontinental resolution (e.g., LA1 vs LA2)
- Somalier uses a separate, smaller set of ~17K coding-region SNPs and does NOT produce continuous ancestry proportions; its role is sample QC and relatedness checking

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `graf_snp_bed` | `null` | GRAF ancestry-informative SNP positions (BED) |
| `grafanc_data` | `null` | GRAF-anc data directory |
| `ancestry_reference_panel` | `null` | Reference panel VCF for supervised classification |
| `ancestry_reference_labels` | `null` | Population labels for reference panel |
| `graf_ancestry_categories` | `EUR,AFR_AM,LA1,LA2,EAS,SAS` | Categories for classification |
| `ancestry_continuous_vars` | `pct_african,...,pct_european` | Continuous proportion variable names |
| `somalier_sites` | `null` | Somalier sites VCF |

**Outputs:** Per-sample `*.graf_snps.vcf.gz` and `*.graf_allele_counts.tsv`; `grafanc_results.txt`; `ancestry_proportions.tsv` (5 continuous components); `ancestry_categories.tsv` (6 discrete groups); Somalier `somalier.samples.tsv`, `somalier.pairs.tsv`, `somalier.html`.

**Software:** bcftools 1.19, GRAF-anc, Somalier 0.2.19, scikit-learn, R

---

### 3.2 Variant Calling and TMB

**Module:** `modules/local/variant_calling/main.nf`

SNV/indel calling from RNA-seq following GATK best practices for RNA-seq, including split N CIGAR reads, base quality score recalibration, hard filtering, and tumor mutational burden estimation.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `GATK_SPLITNCIGARREADS` | Split reads spanning splice junctions into exon-level segments. Reassigns mapping qualities from 255 (STAR default) to 60. Required preprocessing for RNA-seq variant calling. | `process_medium` |
| `GATK_BASERECALIBRATOR` | Model systematic base quality errors using known SNP sites (dbSNP). Generates a recalibration table. Only runs when `known_snps` is provided. | `process_medium` |
| `GATK_APPLYBQSR` | Apply recalibrated base quality scores to the BAM file using the recalibration table. | `process_medium` |
| `GATK_HAPLOTYPECALLER` | Call SNVs and indels via local de novo assembly. Runs in RNA-seq mode with soft-clipped bases excluded. Annotates with dbSNP if provided. | `process_high` |
| `GATK_VARIANTFILTRATION` | Hard-filter variants on: Fisher strand bias (FS > 30), quality by depth (QD < 2), mapping quality (MQ < 40), and depth thresholds. | `process_low` |
| `TMB_ESTIMATION` | Calculate tumor mutational burden from filtered variants per sample. Reports mutations per megabase of callable region. | `process_low` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_base_quality` | `20` | Minimum base quality for variant calling |
| `min_variant_depth` | `10` | Minimum read depth at variant site |
| `min_vaf` | `0.05` | Minimum variant allele fraction |
| `known_snps` / `known_snps_tbi` | `null` | dbSNP VCF and index for BQSR |
| `dbsnp` | `null` | dbSNP VCF for HaplotypeCaller annotation |

**Outputs:** `*.raw.vcf.gz` (raw calls), `*.filtered.vcf.gz` (hard-filtered), `tmb_scores.tsv` and `tmb_summary.tsv` (per-sample TMB).

**Software:** GATK 4.5.0.0, R

---

### 3.3 HLA Typing

**Module:** `modules/local/hla_typing/main.nf`

HLA class I and class II allele typing from RNA-seq BAMs. HLA types feed directly into neoantigen prediction.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `ARCASHLA_EXTRACT` | Extract HLA-mapped reads from BAM into paired FASTQ files targeting the HLA gene region. | `process_medium` |
| `ARCASHLA_GENOTYPE` | Genotype HLA alleles to 4-digit resolution for up to 7 loci: A, B, C (class I) and DPB1, DQB1, DQA1, DRB1 (class II). | `process_medium` |
| `ARCASHLA_MERGE` | Merge per-sample HLA genotype JSONs into a single cohort-level TSV. | `process_low` |
| `OPTITYPE` | HLA class I typing (A, B, C only) using integer linear programming. Often considered a gold standard for class I. Running both tools allows cross-validation of class I calls. | `process_medium` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `hla_tool` | `arcashla` | Tool: `arcashla`, `optitype`, or `both` |
| `arcashla_genes` | `A,B,C,DPB1,DQB1,DQA1,DRB1` | HLA genes to genotype |

**Outputs:** Per-sample `*.genotype.json` (arcasHLA); merged `hla_genotypes.tsv`; OptiType `*_result.tsv` with class I calls.

**Software:** arcasHLA 0.6.0, OptiType 1.3.5

---

## 4. Subworkflow 3: Expression Analysis

**File:** `subworkflows/local/expression_analysis.nf`

Differential expression, co-expression networks, and functional pathway enrichment.

### 4.1 Differential Expression

**Module:** `modules/local/differential_expression/main.nf`

Dual-engine DE analysis using both DESeq2 and limma-voom for concordance validation. Supports all contrast types defined in the contrasts JSON: categorical (two-group), continuous (regression), and categorical_multi (one-vs-reference per level). All contrasts are run within cytomolecular subgroups when specified.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `DESEQ2_DE` | Negative binomial GLM-based differential expression via DESeq2. Adjusts for configurable covariates (batch, sex, age, blast_percentage, tumor_purity, timepoint). Produces per-contrast result tables, volcano plots, MA plots, and heatmaps of top DE genes. Automatically incorporates ancestry proportions when ancestry-related contrasts are present. | `process_high` |
| `LIMMA_VOOM_DE` | Precision-weighted linear modeling on voom-transformed counts via limma. Identical covariate adjustment and contrasts as DESeq2 for direct cross-method comparison. Generates the same suite of per-contrast outputs. | `process_high` |

**Contrast types:**
- **`categorical`**: Two-group comparison (e.g., relapse vs no_relapse) controlling for covariates
- **`continuous`**: Continuous predictor in the DE model (e.g., ancestry proportions); identifies genes scaling linearly with the variable
- **`categorical_multi`**: Auto-generates one contrast per non-reference level (e.g., GRAF categories: AFR_AM vs EUR, LA1 vs EUR, LA2 vs EUR, EAS vs EUR, SAS vs EUR)

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `de_tool` | `both` | Run `deseq2`, `limma`, or `both` |
| `padj_threshold` | `0.05` | Adjusted p-value cutoff |
| `lfc_threshold` | `0.585` | log2 fold-change cutoff (log2(1.5)) |
| `de_covariates` | `batch,sex,age,blast_percentage,tumor_purity,timepoint` | Covariates in the model formula |
| `de_contrasts` | `null` | Path to contrasts JSON (falls back to `assets/default_contrasts.json`) |
| `cytomolecular_subgroups` | `null` | Subgroups for within-group DE (e.g., `ETV6-RUNX1,Hyperdiploid,BCR-ABL1`) |

**Outputs:** Per-contrast `deseq2_results/` and `limma_results/` directories containing result tables (gene, log2FC, pvalue, padj), volcano plots, MA plots, heatmaps, PCA plots; serialized R objects.

**Software:** DESeq2, limma, edgeR (Bioc 3.19)

---

### 4.2 WGCNA Co-expression Networks

**Module:** `modules/local/wgcna/main.nf`

Identifies modules of co-expressed genes and correlates them with sample traits including ancestry proportions, clinical variables, and molecular features.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `WGCNA_ANALYSIS` | Construct a signed weighted co-expression network: (1) auto-detect soft-thresholding power from scale-free topology fit, (2) compute topological overlap matrix (TOM), (3) detect gene modules via dynamic tree cutting, (4) compute module eigengenes, (5) test module-trait correlations (ancestry proportions, relapse, subtype, etc.), (6) identify hub genes per module. | `process_high` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `wgcna_min_module_size` | `30` | Minimum number of genes per module |
| `wgcna_merge_cut_height` | `0.25` | Height cut for merging similar modules |
| `wgcna_soft_power` | `null` | Soft-thresholding power (auto-detected if null) |
| `wgcna_network_type` | `signed` | Network type: `signed` or `unsigned` |

**Outputs:** `module_eigengenes.tsv`, `module_membership.tsv` (gene-to-module with kME scores), `module_trait_cor.tsv` (correlation matrix with p-values), `hub_genes.tsv`; plots: scale-free topology fit, dendrogram with color assignments, module-trait correlation heatmap.

**Software:** WGCNA (Bioc 3.19)

---

### 4.3 Pathway Enrichment

**Module:** `modules/local/pathway_enrichment/main.nf`

Functional enrichment of DE gene lists using three complementary approaches across six MSigDB pathway databases.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `PATHWAY_ENRICHMENT` | For each DE contrast: (1) Over-representation analysis (ORA) on significant gene lists using clusterProfiler, and (2) gene set enrichment analysis (GSEA) on full ranked gene lists using fgsea. Queries all configured MSigDB collections. Generates dotplots, enrichment maps, ridgeplots, and running enrichment score plots. | `process_high` |
| `GSVA_ANALYSIS` | Compute per-sample gene set variation analysis (GSVA) scores for all configured MSigDB collections. Tests for differential pathway activity across sample groups. Produces a sample-by-pathway score matrix and group comparison statistics. | `process_high` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pathway_databases` | `GO_BP,GO_MF,KEGG,REACTOME,HALLMARK,IMMUNESIGDB` | MSigDB collections |
| `msigdb_species` | `Homo sapiens` | Species for msigdbr |
| `gsea_min_size` | `15` | Minimum gene set size for GSEA |
| `gsea_max_size` | `500` | Maximum gene set size |
| `run_gsva` | `true` | Enable GSVA analysis |

**Outputs:** `pathway_results/ora/` and `pathway_results/gsea/` per contrast and database; GSEA NES, p-values, and leading edge genes; `gsva_scores.tsv` (sample-by-pathway matrix); pathway and GSVA plots.

**Software:** clusterProfiler, fgsea, GSVA, msigdbr (Bioc 3.19 / CRAN)

---

## 5. Subworkflow 4: Structural Variants

**File:** `subworkflows/local/structural_variants.nf`

Differential and aberrant splicing, gene fusion detection, and expression-based CNV inference.

### 5.1 Differential Splicing

**Module:** `modules/local/splicing/main.nf`

Four complementary splicing analysis approaches, each capturing different aspects of alternative splicing:

**rMATS** -- Event-level statistical testing:

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `RMATS` | Detect differential splicing events between two groups via a likelihood-ratio test. Quantifies five event types: skipped exon (SE), alternative 3' splice site (A3SS), alternative 5' splice site (A5SS), mutually exclusive exons (MXE), and retained intron (RI). Reports inclusion level differences (delta PSI) and significance. | `process_high` |

**Leafcutter** -- Intron cluster-level analysis:

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `LEAFCUTTER_JUNCTIONS` | Extract splice junctions from BAM using regtools. | `process_low` |
| `LEAFCUTTER_CLUSTER` | Cluster introns across samples into groups sharing a splice site. | `process_medium` |
| `LEAFCUTTER_DIFF_SPLICING` | Dirichlet-multinomial GLM test for differential intron usage between groups. Annotation-free approach that does not require predefined event types. | `process_medium` |

**SplAdder** -- Graph-based event detection:

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `SPLADDER_BUILD` | Construct augmented splicing graphs from BAMs + GTF annotation. Detect and quantify alternative splicing events with configurable confidence levels. | `process_high` |
| `SPLADDER_TEST` | Negative binomial GLM differential testing of SplAdder-detected events between sample groups. | `process_medium` |

**Bisbee** -- Beta-binomial modeling with protein effects:

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `BISBEE_PREP` | Extract inclusion/exclusion junction counts from SplAdder HDF5 output files for Bisbee input. | `process_medium` |
| `BISBEE_DIFF` | Beta-binomial differential splicing test that accounts for biological overdispersion. More robust than binomial tests for small sample sizes. | `process_medium` |
| `BISBEE_PROT` | Predict protein-level effects of differential splice events: nonsense-mediated decay (NMD), frameshifts, protein domain disruptions. | `process_medium` |
| `BISBEE_OUTLIER` | Per-sample splicing outlier detection using a beta-binomial model. Identifies individual samples with aberrant splicing at specific junctions (e.g., driver splice mutations). | `process_medium` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `splicing_tool` | `all` | Tools to run: `rmats`, `leafcutter`, `spladder`, `bisbee`, `both`, or `all` |
| `rmats_read_length` | `150` | Read length for rMATS |
| `rmats_novel_ss` | `true` | Allow novel splice sites in rMATS |
| `leafcutter_min_coverage` | `20` | Minimum junction read support for Leafcutter clustering |
| `spladder_confidence` | `3` | SplAdder confidence level (1-3; higher = stricter) |
| `bisbee_outlier_fdr` | `0.05` | FDR threshold for Bisbee outlier detection |

**Outputs:** rMATS per-event-type result tables; Leafcutter cluster significance and effect sizes; SplAdder splicing graphs and test results; Bisbee differential results, protein effect predictions (NMD targets, frameshifts, domain disruptions), and per-sample outlier calls.

**Software:** rMATS 4.3.0, Leafcutter 0.2.9, regtools, SplAdder 3.0.4, Bisbee

---

### 5.2 Gene Fusion Detection

**Module:** `modules/local/fusions/main.nf`

Dual-caller fusion detection with multi-caller merging and visualization.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `FUSIONCATCHER` | Detect fusions via split-read and paired-end mapping. Converts BAM to FASTQ internally. Uses built-in database of known oncogenic fusions for annotation. | `process_high` |
| `ARRIBA` | Detect fusions from STAR-aligned BAMs using split-read and discordant-pair evidence. Supports blacklist filtering of recurrent artifacts and known-fusion boosting. | `process_high` |
| `ARRIBA_VISUALIZATION` | Generate publication-quality PDF plots of fusion breakpoints with protein domain annotations and exon structure. | `process_medium` |
| `MERGE_FUSIONS` | Merge and prioritize calls from both callers. Fusions detected by multiple tools are flagged as high-confidence. Reports per-sample fusion counts and cross-sample recurrence. | `process_low` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fusion_tool` | `both` | Run `fusioncatcher`, `arriba`, or `both` |
| `min_fusion_reads` | `3` | Minimum supporting reads to report a fusion |
| `fusioncatcher_data` | `null` | Path to FusionCatcher reference database |
| `arriba_blacklist` | `null` | Arriba blacklist file |
| `arriba_known_fusions` | `null` | Known fusions for Arriba annotation |

**Outputs:** FusionCatcher `final-list_candidate-fusion-genes.txt`; Arriba `*.arriba.fusions.tsv` and `*.arriba.pdf` visualizations; `merged_fusions.tsv`; `high_confidence_fusions.tsv`.

**Software:** FusionCatcher 1.33, Arriba 2.4.0, Samtools 1.19

---

### 5.3 CNV Inference

**Module:** `modules/local/cnv_inference/main.nf`

Infer large-scale copy number variations from RNA-seq expression data.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `INFERCNV` | Compare tumor gene expression profiles across chromosomal positions against a normal reference group. Applies noise filtering and HMM-based state prediction to infer chromosomal gains, losses, and focal events. Generates genome-wide expression heatmaps with inferred CNV states. | `process_very_high` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cnv_tool` | `infercnv` | CNV inference tool |
| `cnv_reference_group` | `null` | Normal reference group label in annotations file |
| `cnv_gene_order_file` | `null` | Gene position ordering file (gene, chr, start, end) |

**Notes:**
- Computationally intensive; labeled `process_very_high` (16 CPUs, 128 GB)
- Requires a reference group of normal/non-tumor samples
- Particularly useful for detecting B-ALL-associated aneuploidies (hyperdiploidy, hypodiploidy)

**Outputs:** InferCNV results directory, `infercnv.png` heatmap, `infercnv.observations.txt` per-gene CNV scores.

**Software:** inferCNV (Bioc 3.19)

---

## 6. Subworkflow 5: Immunogenomics

**File:** `subworkflows/local/immunogenomics.nf`

Immune cell deconvolution, neoantigen prediction from three RNA-seq-specific sources, and TCR/BCR repertoire profiling.

### 6.1 Immune Deconvolution

**Module:** `modules/local/immune_analysis/main.nf`

Estimates immune cell composition and tumor purity from expression data using multiple deconvolution methods.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `IMMUNE_DECONVOLUTION` | Run up to six deconvolution methods via the immunedeconv R framework: CIBERSORTx (LM22 signature matrix), xCell (64 cell types via ssGSEA), MCP-counter (10 immune/stromal populations), EPIC (6 major immune + other), TIMER (6 tumor-infiltrating immune cell types), and quantiseq. Optionally corrects cell fractions for tumor purity. | `process_high` |
| `ESTIMATE_SCORES` | Compute stromal score, immune score, ESTIMATE score, and inferred tumor purity using the ESTIMATE algorithm. Provides an independent computational estimate of tumor purity from expression data. | `process_medium` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `immune_deconv_methods` | `cibersortx,xcell,mcpcounter,estimate,epic,timer` | Methods to run |
| `cibersortx_token` | `null` | CIBERSORTx API token (required for CIBERSORTx) |
| `cibersortx_sigmatrix` | `LM22` | Signature matrix for CIBERSORTx |
| `estimate_platform` | `illumina` | Platform type for ESTIMATE |
| `correct_tumor_purity` | `true` | Adjust immune fractions for tumor purity |

**Outputs:** `deconvolution_all.tsv` (all methods combined), `cell_fractions.tsv`, `estimate_scores.tsv`, `estimate_purity.tsv`; box plots, heatmaps, correlation plots.

**Software:** immunedeconv, ESTIMATE, CIBERSORTx, xCell, MCP-counter, EPIC, TIMER (Bioc 3.19)

---

### 6.2 Neoantigen Prediction

**Module:** `modules/local/neoantigen/main.nf`

RNA-seq-specific neoantigen prediction from three mutation sources -- no matched WGS/WES required. Predicts MHC class I and class II binding using patient-specific HLA types from the HLA typing module.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `PVACSEQ` | Predict neoantigens from somatic SNVs/indels. Evaluates MHC class I (8-11mer) and class II (15mer) binding via MHCflurry, NetMHCpan, and NetMHCIIpan. Filters on binding affinity and percentile rank thresholds. | `process_high` |
| `NEOFUSE` | Predict fusion-derived neoantigens from Arriba fusion calls using patient HLA alleles. Evaluates novel peptide sequences spanning fusion breakpoints. | `process_high` |
| `SNAF_SPLICING_NEOANTIGENS` | Identify neoantigens from alternative splicing junctions. Predicts both T-cell (MHC-restricted) and B-cell (surface-accessible) epitopes from novel splice-derived peptides. Provides orthogonal neoantigen detection to pVACseq and NeoFuse. | `process_high` |
| `MERGE_NEOANTIGENS` | Merge and deduplicate candidates from all three sources. Compute per-sample neoantigen burden stratified by source (SNV, fusion, splicing). | `process_low` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `neoantigen_tool` | `all` | Tools: `pvacseq`, `neofuse`, `snaf`, `both`, or `all` |
| `pvac_algorithms` | `MHCflurry,NetMHCpan,NetMHCIIpan` | Binding prediction algorithms |
| `binding_threshold` | `500` | Binding affinity threshold (nM) |
| `percentile_threshold` | `2.0` | Percentile rank cutoff |
| `snaf_db` | `null` | Path to SNAF reference database |

**Outputs:** Per-sample pVACseq, NeoFuse, and SNAF results; `merged_neoantigens.tsv`; `neoantigen_burden.tsv`.

**Software:** pVACtools (pVACseq) 4.2.0, NeoFuse 1.0, SNAF

---

### 6.3 TCR/BCR Repertoire

**Module:** `modules/local/tcr_repertoire/main.nf`

Reconstruct adaptive immune receptor repertoires from bulk RNA-seq reads.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `TRUST4` | Extract and assemble TCR/BCR CDR3 sequences from RNA-seq BAM reads mapping to immune receptor loci. Reconstructs both TCR (alpha/beta, gamma/delta) and BCR (heavy/light) chains. Reports clonotype frequency, V/D/J gene usage, and CDR3 amino acid sequences. | `process_high` |
| `MERGE_TCR_REPORTS` | Merge per-sample TRUST4 reports into cohort-level summaries. Compute diversity indices (Shannon entropy, Simpson index, clonality) and track shared clonotypes across samples. | `process_low` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tcr_tool` | `trust4` | TCR/BCR repertoire tool |

**Notes:**
- Sensitivity depends on immune cell infiltration; bone marrow samples (B-ALL) typically have higher yields
- Diversity metrics enable comparison of immune repertoire complexity between groups

**Outputs:** Per-sample clonotype reports; `tcr_repertoire_summary.tsv`; `tcr_diversity_metrics.tsv` (Shannon, Simpson, clonality); `tcr_clonotype_tracking.tsv`.

**Software:** TRUST4 1.1.0

---

## 7. Subworkflow 6: Clinical Analysis

**File:** `subworkflows/local/clinical_analysis.nf`

Sensitivity assessment, pharmacological target identification, and molecular subtype classification.

### 7.1 Sensitivity Analysis

**Module:** `modules/local/sensitivity/main.nf`

Assesses the robustness of DE results to analytical choices and model specifications.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `SENSITIVITY_ANALYSIS` | Four systematic assessments: (1) **Timepoint-stratified DE** -- re-run DE within diagnostic-only and relapse-only samples; (2) **Ancestry threshold sensitivity** -- vary ancestry proportion cutoffs and measure DE result stability; (3) **Covariate leave-one-out** -- remove each covariate individually and measure impact on DE gene lists; (4) **Model comparison** -- DESeq2 vs limma concordance under different model specifications. Uses bootstrap resampling (default 1000 iterations) for stability estimates. | `process_high` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sensitivity_variables` | `timepoint,ancestry_proportion` | Variables to assess |
| `bootstrap_iterations` | `1000` | Bootstrap iterations for stability |

**Outputs:** `timepoint_analysis.tsv`, `ancestry_sensitivity.tsv`, `covariate_impact.tsv`, `model_comparison.tsv`; Upset plots, Jaccard similarity heatmaps, forest plots.

**Software:** R (DESeq2, limma) (Bioc 3.19)

---

### 7.2 Pharmacogenomics

**Module:** `modules/local/pharmacogenomics/main.nf`

Predicts drug sensitivity from expression profiles and identifies druggable targets.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `PHARMACOGENOMICS` | Four-part analysis: (1) **Drug sensitivity prediction** -- per-sample IC50/AUC scores using oncoPredict trained on GDSC/CCLE cell line pharmacogenomics data; (2) **Druggable target identification** -- query DGIdb for drug-gene interactions from DE gene lists; (3) **CMap connectivity** (optional) -- query Connectivity Map signatures for drug repurposing; (4) **Group comparisons** -- compare predicted drug sensitivity between ancestry and clinical groups. | `process_high` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `drug_response_db` | `GDSC` | Training database: `GDSC`, `CCLE`, `PRISM`, or `all` |
| `dgidb_interactions` | `true` | Query DGIdb for druggable targets |
| `cmap_signatures` | `null` | Path to CMap signature database |

**Outputs:** `drug_sensitivity_scores.tsv`, `druggable_targets.tsv`, `dgidb_interactions.tsv`, `cmap_connections.tsv` (optional), `group_comparisons/`; drug sensitivity heatmaps, waterfall plots.

**Software:** oncoPredict (Bioc 3.19), DGIdb

---

### 7.3 Molecular Subtyping

**Module:** `modules/local/molecular_subtyping/main.nf`

Expression-based classification into disease subtypes.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `MOLECULAR_SUBTYPING` | Classify each sample into a molecular subtype using gene expression signatures and consensus scoring across multiple gene set methods (e.g., GSVA, ssGSEA, z-score). Outputs subtype assignments with confidence scores. | `process_medium` |

**Default B-ALL subtypes** (12+ signatures):

| Subtype | Cytogenetic Feature |
|---------|-------------------|
| ETV6-RUNX1 | t(12;21) |
| BCR-ABL1 / Ph-like | t(9;22) or Ph-like kinase activation |
| KMT2A-rearranged | MLL rearrangements |
| High hyperdiploidy | 51-67 chromosomes |
| Low hypodiploidy / near-haploid | <44 chromosomes |
| TCF3-PBX1 | t(1;19) |
| iAMP21 | Intrachromosomal amplification of chromosome 21 |
| DUX4-rearranged | DUX4 insertions into IGH |
| MEF2D-rearranged | MEF2D fusions |
| ZNF384-rearranged | ZNF384 fusions |
| NUTM1-rearranged | NUTM1 fusions |
| PAX5alt / PAX5 P80R | PAX5 alterations |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `subtyping_method` | `consensus` | Classification method |

**Outputs:** `subtypes.tsv` (per-sample assignments), `classifier_scores.tsv` (per-subtype scores), heatmaps and UMAP/PCA plots by subtype.

**Software:** R (Bioc 3.19)

---

## 8. Validation

**Module:** `modules/local/validation/main.nf`

Cross-pipeline quality control checks and artifact detection.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `VALIDATE_DE_CONCORDANCE` | Compare DESeq2 and limma-voom results per contrast: Spearman rank correlation of log2FC, Jaccard index of significant gene overlap, log2FC agreement scatter, and discordant gene flagging. | `process_medium` |
| `VALIDATE_SEX_CHECK` | Verify reported sample sex against XIST expression and Y-chromosome gene expression (RPS4Y1, EIF1AY, DDX3Y, KDM5D). Flag mismatches indicating potential sample swaps. | `process_low` |
| `VALIDATE_EXPRESSION_OUTLIERS` | Detect outlier samples using: PCA distance from centroid (configurable SD threshold), pairwise sample correlation (minimum threshold), and expression distribution metrics. | `process_medium` |
| `VALIDATE_GENOMIC_INFLATION` | Compute genomic inflation factor (lambda) from DE p-value distributions per contrast. Lambda > 1.1 indicates systematic bias. Generate QQ plots for visual assessment. | `process_low` |
| `COMPILE_VALIDATION_REPORT` | Aggregate all validation findings into a unified markdown report with pass/fail summary table. | `process_low` |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `outlier_sd_threshold` | `3` | SD threshold for PCA-based outlier detection |
| `min_sample_correlation` | `0.8` | Minimum pairwise correlation to pass QC |

**Outputs:** `concordance_summary.tsv`, `sex_check/` results, `expression_outliers/` list, lambda values and QQ plots per contrast, compiled `validation_report.md` and `validation_summary.tsv`.

**Software:** R (Bioc 3.19)

---

## 9. Visualization and Reporting

**Module:** `modules/local/visualization/main.nf`

Publication-quality figure generation and pipeline analysis logging.

| Process | Description | Resource Label |
|---------|-------------|----------------|
| `VISUALIZATION_SUMMARY` | Create integrated summary figures spanning DE (volcano, MA, heatmaps), pathway (dotplots, running enrichment), immune (composition bars, score heatmaps), ancestry (PCA, stacked bars), WGCNA (dendrogram, module-trait heatmap), and pharmacogenomics (drug sensitivity heatmap). PNG + PDF formats. Generates an HTML figure index for browsing. | `process_medium` |
| `PIPELINE_LOGGER` | Record a structured log entry for each pipeline step with timestamp, description, input files, and key findings. Each module appends to a shared `findings.md`. | `process_low` |
| `COMPILE_PIPELINE_LOG` | Aggregate all log entries and module findings into a comprehensive `pipeline_analysis_log.md` with pipeline parameters, step summaries, and flagged results. | `process_low` |

**Figure output structure:**
```
figures/
├── de/           # Volcano, MA, heatmaps per contrast
├── pathway/      # Dotplots, enrichment scores, GSVA heatmap
├── immune/       # Composition bars, score heatmaps, purity plots
├── ancestry/     # PCA scatters, proportion stacked bars
├── wgcna/        # Dendrogram, module-trait heatmap, eigengene plots
├── pharma/       # Drug sensitivity heatmap, waterfall plots
├── overview/     # Cross-module summary figures
└── figure_index.html
```

**Pipeline logging captures:**
- QC flags (samples failing thresholds, strandedness mismatches)
- Ancestry assignments and SNP coverage statistics
- Genome build detections and liftover outcomes
- DE summary (significant gene counts, cross-method concordance)
- Splicing outliers (Bisbee per-sample flags)
- Fusion calls (high-confidence multi-caller fusions)
- Neoantigen counts by source (SNV, fusion, splicing)
- Validation warnings (sex mismatches, PCA outliers, lambda > 1.1)

**Software:** R (ggplot2, ComplexHeatmap), Python

---

## 10. Complete Software Reference

All tools used across the pipeline, with versions and container sources defined in `conf/containers.config`.

| Tool | Version | Category | Used By |
|------|---------|----------|---------|
| Samtools | 1.19 | Alignment utilities | QC, genome utils, fusions |
| RSeQC | 5.0.3 | RNA-seq QC | QC |
| Picard | 3.1.1 | Alignment QC | QC |
| MultiQC | 1.21 | QC aggregation | QC |
| CrossMap | -- | Coordinate liftover | Genome utils |
| Subread (featureCounts) | 2.0.6 | Read counting | Counting |
| DESeq2 | Bioc 3.19 | Differential expression | Counting, DE, sensitivity |
| limma | Bioc 3.19 | Differential expression | DE, sensitivity |
| edgeR | Bioc 3.19 | Normalization | Counting, DE |
| sva (ComBat-seq) | Bioc 3.19 | Batch correction | Counting |
| bcftools | 1.19 | Variant utilities | Ancestry |
| GRAF-anc | -- | Ancestry inference | Ancestry |
| Somalier | 0.2.19 | Sample QC/relatedness | Ancestry |
| scikit-learn | -- | Classification | Ancestry |
| GATK | 4.5.0.0 | Variant calling | Variant calling |
| arcasHLA | 0.6.0 | HLA typing | HLA typing |
| OptiType | 1.3.5 | HLA typing | HLA typing |
| rMATS | 4.3.0 | Differential splicing | Splicing |
| Leafcutter | 0.2.9 | Differential splicing | Splicing |
| regtools | -- | Junction extraction | Splicing |
| SplAdder | 3.0.4 | Splicing graphs | Splicing |
| Bisbee | -- | Splicing + protein effects | Splicing |
| FusionCatcher | 1.33 | Fusion detection | Fusions |
| Arriba | 2.4.0 | Fusion detection | Fusions |
| inferCNV | Bioc 3.19 | CNV inference | CNV |
| WGCNA | Bioc 3.19 | Co-expression networks | WGCNA |
| clusterProfiler | Bioc 3.19 | Pathway enrichment (ORA) | Pathway |
| fgsea | Bioc 3.19 | Pathway enrichment (GSEA) | Pathway |
| GSVA | Bioc 3.19 | Pathway activity | Pathway |
| msigdbr | CRAN | MSigDB gene sets | Pathway |
| immunedeconv | Bioc 3.19 | Immune deconvolution | Immune |
| ESTIMATE | Bioc 3.19 | Tumor purity | Immune |
| CIBERSORTx | -- | Immune deconvolution | Immune |
| xCell | -- | Immune deconvolution | Immune |
| MCP-counter | -- | Immune deconvolution | Immune |
| EPIC | -- | Immune deconvolution | Immune |
| TIMER | -- | Immune deconvolution | Immune |
| pVACtools (pVACseq) | 4.2.0 | Neoantigen prediction | Neoantigen |
| NeoFuse | 1.0 | Fusion neoantigens | Neoantigen |
| SNAF | -- | Splicing neoantigens | Neoantigen |
| TRUST4 | 1.1.0 | TCR/BCR repertoire | TCR repertoire |
| oncoPredict | Bioc 3.19 | Drug sensitivity | Pharmacogenomics |
| DGIdb | -- | Drug-gene interactions | Pharmacogenomics |

---

*Generated for the Cancer Bulk RNA-Seq Analysis Pipeline with Ancestry Integration v1.0.0*
