# Cancer Bulk RNA-Seq Analysis Pipeline with Ancestry Integration

A comprehensive Nextflow DSL2 pipeline for cancer bulk RNA-seq analysis with integrated genetic ancestry inference, immunogenomics, and pharmacogenomics. Default covariates and molecular subtyping signatures are tuned for pediatric B-cell acute lymphoblastic leukemia (B-ALL), but the pipeline is fully configurable for any cancer type.

**Version:** 1.0.0 | **Nextflow:** >=23.04.0 | **DSL:** 2

---

## Table of Contents

1. [Overview](#1-overview)
2. [Pipeline Diagram](#2-pipeline-diagram)
3. [Quick Start](#3-quick-start)
4. [Modules](#4-modules)
5. [Samplesheet Format](#5-samplesheet-format)
6. [Contrast Definitions](#6-contrast-definitions)
7. [Covariates](#7-covariates)
8. [Differential Variables](#8-differential-variables)
9. [Parameters](#9-parameters)
10. [Output Directory Structure](#10-output-directory-structure)
11. [Configurability for Other Cancers](#11-configurability-for-other-cancers)
12. [Profiles](#12-profiles)
13. [Pipeline Logging](#13-pipeline-logging)
14. [Visualizations Produced](#14-visualizations-produced)
15. [Credits and References](#15-credits-and-references)

---

## 1. Overview

| Feature | Detail |
|---|---|
| **Modules** | 20+ analysis modules organized into 6 subworkflows |
| **Input** | STAR-aligned, sorted, duplicate-marked BAM/BAI files |
| **Genome build** | Auto-detects hg19 vs hg38 from BAM header (chr1 length); CrossMap liftover when needed |
| **DE approach** | Dual engine (DESeq2 + limma-voom) with concordance validation |
| **Ancestry** | GRAF-pop 282K ancestry-informative SNPs; continuous and categorical ancestry as covariates |
| **Neoantigen** | RNA-seq only pipeline (pVACseq, NeoFuse, SNAF) -- no matched WGS required |
| **Containers** | Every process containerized; Docker, Singularity, and Conda supported |
| **Schedulers** | SLURM, LSF, local, and custom profiles |
| **Default disease** | Pediatric B-ALL (12+ cytomolecular subtype signatures) |

### Key Capabilities

- Quality control with samtools, RSeQC, Picard, and MultiQC
- Genetic ancestry inference from RNA-seq reads (continental and subcontinental)
- Read counting, normalization (VST/rlog/TMM/TPM), and batch correction (ComBat-seq)
- Variant calling from RNA-seq (GATK best practices) with TMB estimation
- Differential expression with ancestry and clinical covariates
- Differential splicing (rMATS, Leafcutter, SplAdder, Bisbee)
- Gene fusion detection (Arriba, FusionCatcher)
- CNV inference from expression (InferCNV)
- Weighted gene co-expression networks (WGCNA)
- Pathway enrichment (ORA, GSEA, GSVA) across MSigDB collections
- Immune deconvolution (xCell, MCP-counter, EPIC, CIBERSORTx, TIMER) corrected for tumor purity
- HLA typing (arcasHLA, OptiType) and neoantigen prediction
- TCR/BCR repertoire profiling (TRUST4)
- Molecular subtyping via curated gene signatures
- Pharmacogenomics (oncoPredict/GDSC, DGIdb)
- Sensitivity and validation analyses

---

## 2. Pipeline Diagram

```
                             STAR-aligned BAM/BAI files
                                       |
                                       v
                          +------------------------+
                          |  Genome Build Detection |
                          |  (hg19 auto-liftover)  |
                          +------------------------+
                                       |
                    +------------------+------------------+
                    |                                     |
                    v                                     v
          +------------------+                  +------------------+
          |   QC & Metrics   |                  |   featureCounts  |
          | samtools, RSeQC, |                  |   + Normalize    |
          | Picard, MultiQC  |                  |   + Batch Corr.  |
          +------------------+                  +------------------+
                    |                                     |
                    +------------------+------------------+
                                       |
          +----------------------------+----------------------------+
          |              |             |             |              |
          v              v             v             v              v
  +-----------+  +-----------+  +-----------+  +-----------+  +-----------+
  | GENOMICS  |  |EXPRESSION |  |STRUCTURAL |  | IMMUNO-   |  | CLINICAL  |
  |           |  | ANALYSIS  |  | VARIANTS  |  | GENOMICS  |  | ANALYSIS  |
  +-----------+  +-----------+  +-----------+  +-----------+  +-----------+
  | Ancestry  |  | DESeq2    |  | rMATS     |  | Immune    |  | Sensitiv. |
  | GRAF-anc  |  | limma-voom|  | Leafcutter|  |  deconv   |  | Pharma-   |
  | Somalier  |  | WGCNA     |  | SplAdder  |  | ESTIMATE  |  |  genomics |
  | GATK HC   |  | Pathway   |  | Bisbee    |  | pVACseq   |  | Molecular |
  | HLA typing|  |  enrichm. |  | Arriba    |  | NeoFuse   |  |  subtype  |
  | TMB est.  |  | GSVA      |  | FusionC.  |  | SNAF      |  |           |
  |           |  |           |  | InferCNV  |  | TRUST4    |  |           |
  +-----------+  +-----------+  +-----------+  +-----------+  +-----------+
          |              |             |             |              |
          +----------------------------+----------------------------+
                                       |
                                       v
                          +------------------------+
                          |  Validation & Reporting |
                          |  DE concordance, QQ,    |
                          |  sex check, figures,    |
                          |  findings.md, MultiQC   |
                          +------------------------+
```

### Subworkflow Organization

| # | Subworkflow | File | Description |
|---|---|---|---|
| 1 | `PREPROCESSING` | `subworkflows/local/preprocessing.nf` | QC, genome build detection, liftover, featureCounts, normalization, batch correction |
| 2 | `GENOMICS` | `subworkflows/local/genomics.nf` | Ancestry inference, GATK variant calling, HLA typing, TMB estimation |
| 3 | `EXPRESSION_ANALYSIS` | `subworkflows/local/expression_analysis.nf` | DESeq2/limma-voom DE, WGCNA, pathway enrichment, GSVA |
| 4 | `STRUCTURAL_VARIANTS` | `subworkflows/local/structural_variants.nf` | Differential splicing, gene fusions, CNV inference |
| 5 | `IMMUNOGENOMICS` | `subworkflows/local/immunogenomics.nf` | Immune deconvolution, ESTIMATE, neoantigens, TCR/BCR |
| 6 | `CLINICAL_ANALYSIS` | `subworkflows/local/clinical_analysis.nf` | Sensitivity analysis, pharmacogenomics, molecular subtyping |

---

## 3. Quick Start

### Minimal run

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta /refs/GRCh38.fa \
    --fasta_fai /refs/GRCh38.fa.fai \
    --gtf /refs/gencode.v44.annotation.gtf \
    --outdir ./results \
    -profile docker
```

### Full run with ancestry and neoantigen prediction

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta /refs/GRCh38.fa \
    --fasta_fai /refs/GRCh38.fa.fai \
    --gtf /refs/gencode.v44.annotation.gtf \
    --known_snps /refs/dbsnp_146.hg38.vcf.gz \
    --known_snps_tbi /refs/dbsnp_146.hg38.vcf.gz.tbi \
    --graf_snp_bed /refs/graf_snps.282K.bed \
    --grafanc_data /refs/grafanc_data/ \
    --somalier_sites /refs/somalier_sites.hg38.vcf.gz \
    --de_contrasts /path/to/contrasts.json \
    --de_covariates 'batch,sex,age,blast_percentage,tumor_purity,timepoint' \
    --outdir ./results \
    -profile singularity,slurm \
    -resume
```

### Test run (CI/CD)

```bash
nextflow run main.nf -profile test,docker
```

---

## 4. Modules

### 4.1 QC and BAM Metrics

| Process | Tool | Description |
|---|---|---|
| `SAMTOOLS_FLAGSTAT` | [samtools](http://www.htslib.org/) 1.19 | Alignment flag statistics |
| `SAMTOOLS_IDXSTATS` | [samtools](http://www.htslib.org/) 1.19 | Per-chromosome read counts |
| `SAMTOOLS_STATS` | [samtools](http://www.htslib.org/) 1.19 | Comprehensive alignment statistics |
| `RSEQC_BAMSTAT` | [RSeQC](https://rseqc.sourceforge.net/) 5.0.3 | RNA-seq BAM quality metrics |
| `RSEQC_READDISTRIBUTION` | [RSeQC](https://rseqc.sourceforge.net/) 5.0.3 | Read distribution across genomic features |
| `RSEQC_INFEREXPERIMENT` | [RSeQC](https://rseqc.sourceforge.net/) 5.0.3 | Strandedness inference |
| `PICARD_COLLECTRNASEQMETRICS` | [Picard](https://broadinstitute.github.io/picard/) 3.1.1 | RNA-seq specific QC metrics |
| `MULTIQC` | [MultiQC](https://multiqc.info/) 1.21 | Aggregate QC report |

**Outputs:** `results/qc/samtools/`, `results/qc/rseqc/`, `results/qc/picard/`, `results/qc/multiqc/`

---

### 4.2 Genome Utilities

| Process | Tool | Description |
|---|---|---|
| `DETECT_GENOME_BUILD` | samtools | Auto-detect hg19 vs hg38 from BAM header (chr1 contig length) |
| `CROSSMAP_BAM` | [CrossMap](https://crossmap.sourceforge.net/) | Liftover hg19 BAMs to hg38 coordinates |

**Outputs:** Lifted-over BAMs fed into downstream modules (intermediate, not published)

---

### 4.3 Counting and Normalization

| Process | Tool | Description |
|---|---|---|
| `SUBREAD_FEATURECOUNTS` | [Subread featureCounts](https://subread.sourceforge.net/) 2.0.6 | Gene-level read counting |
| `MERGE_COUNTS` | R (custom) | Merge per-sample counts into a matrix |
| `NORMALIZE_COUNTS` | R / [DESeq2](https://bioconductor.org/packages/DESeq2/) | VST, rlog, TMM, or TPM normalization |
| `BATCH_CORRECTION` | R / [sva (ComBat-seq)](https://bioconductor.org/packages/sva/) | Remove batch effects from raw counts |

**Parameters:** `fc_strandedness` (0/1/2), `normalization_method` (vst/rlog/tmm/tpm), `batch_correction` (combat_seq/sva/none)

**Outputs:** `results/counting/featurecounts/`, `results/counting/normalized/`, `results/counting/batch_corrected/`

---

### 4.4 Genetic Ancestry Inference

| Process | Tool | Description |
|---|---|---|
| `EXTRACT_GRAF_SNPS` | [bcftools](https://samtools.github.io/bcftools/) 1.19 | Extract allele counts at 282K GRAF ancestry-informative SNP positions |
| `ANCESTRY_INFERENCE` | Python (custom; [GRAF-pop](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/Software.cgi) method) | Supervised classification into continental/subcontinental ancestry categories; outputs continuous proportions |
| `SOMALIER_EXTRACT` | [Somalier](https://github.com/brentp/somalier) 0.2.19 | Extract genotypes at informative sites for QC |
| `SOMALIER_ANCESTRY` | [Somalier](https://github.com/brentp/somalier) 0.2.19 | Sample relatedness and ancestry QC |

**Coverage notes:**
- RNA-seq typically covers 20K-50K of the 282K GRAF SNPs (expressed-region SNPs only)
- 10K+ SNPs: sufficient for continental-level ancestry (EUR/AFR/EAS/SAS/AMR)
- 20K+ SNPs: recommended for subcontinental resolution (e.g., LA1 vs LA2 distinction)

**GRAF ancestry categories:** EUR, AFR_AM, LA1, LA2, EAS, SAS

**Outputs:** `results/ancestry/graf_snps/`, `results/ancestry/proportions/`, `results/ancestry/somalier/`

---

### 4.5 Variant Calling and TMB

| Process | Tool | Description |
|---|---|---|
| `GATK_SPLITNCIGARREADS` | [GATK](https://gatk.broadinstitute.org/) 4.5.0 | Split reads spanning splice junctions |
| `GATK_BASERECALIBRATOR` | [GATK](https://gatk.broadinstitute.org/) 4.5.0 | Base quality score recalibration table |
| `GATK_APPLYBQSR` | [GATK](https://gatk.broadinstitute.org/) 4.5.0 | Apply recalibration to BAM |
| `GATK_HAPLOTYPECALLER` | [GATK](https://gatk.broadinstitute.org/) 4.5.0 | Call variants in RNA-seq mode |
| `GATK_VARIANTFILTRATION` | [GATK](https://gatk.broadinstitute.org/) 4.5.0 | Hard-filter RNA-seq variants |
| `TMB_ESTIMATION` | R (custom) | Tumor mutational burden from filtered variants |

**Outputs:** `results/variant_calling/raw/`, `results/variant_calling/filtered/`, `results/variant_calling/tmb/`

---

### 4.6 Differential Expression

| Process | Tool | Description |
|---|---|---|
| `DESEQ2_DE` | R / [DESeq2](https://bioconductor.org/packages/DESeq2/) | Negative binomial GLM differential expression |
| `LIMMA_VOOM_DE` | R / [limma](https://bioconductor.org/packages/limma/) | Precision-weighted linear model DE |

Both tools are run by default (`de_tool = 'both'`) for concordance validation. Contrasts are defined via a JSON file (see [Section 6](#6-contrast-definitions)). Covariates are fully configurable (see [Section 7](#7-covariates)).

**Outputs:** `results/differential_expression/deseq2/`, `results/differential_expression/limma/`

---

### 4.7 Differential Splicing

| Process | Tool | Description |
|---|---|---|
| `RMATS` | [rMATS](http://rnaseq-mats.sourceforge.net/) 4.3.0 | Detect differential splicing events (SE, MXE, A3SS, A5SS, RI) |
| `LEAFCUTTER_JUNCTIONS` | [Leafcutter](https://davidaknowles.github.io/leafcutter/) 0.2.9 | Extract splice junctions from BAMs |
| `LEAFCUTTER_CLUSTER` | [Leafcutter](https://davidaknowles.github.io/leafcutter/) 0.2.9 | Cluster introns across samples |
| `LEAFCUTTER_DIFF_SPLICING` | [Leafcutter](https://davidaknowles.github.io/leafcutter/) 0.2.9 | Differential intron usage analysis |
| `SPLADDER_BUILD` | [SplAdder](https://github.com/ratschlab/spladder) 3.0.4 | Build augmented splicing graphs |
| `SPLADDER_TEST` | [SplAdder](https://github.com/ratschlab/spladder) 3.0.4 | Statistical testing of alternative splicing |
| `BISBEE_PREP` | Python (custom) | Prepare SplAdder output for Bisbee |
| `BISBEE_DIFF` | R / [Bisbee](https://github.com/tgen/bisbee) | Beta-binomial differential splicing |
| `BISBEE_PROT` | Python (custom) | Predict protein-level effects of splicing events |
| `BISBEE_OUTLIER` | R / [Bisbee](https://github.com/tgen/bisbee) | Outlier splicing event detection (per-sample) |

**Outputs:** `results/splicing/rmats/`, `results/splicing/leafcutter/`, `results/splicing/spladder/`, `results/splicing/bisbee/`

---

### 4.8 Gene Fusion Detection

| Process | Tool | Description |
|---|---|---|
| `FUSIONCATCHER` | [FusionCatcher](https://github.com/ndaniel/fusioncatcher) 1.33 | Fusion gene detection using multiple alignment methods |
| `ARRIBA` | [Arriba](https://arriba.readthedocs.io/) 2.4.0 | Fast and accurate fusion detection from STAR alignments |
| `ARRIBA_VISUALIZATION` | [Arriba](https://arriba.readthedocs.io/) 2.4.0 | Publication-quality fusion visualizations |
| `MERGE_FUSIONS` | Python (custom) | Merge and deduplicate fusions from both callers |

**Outputs:** `results/fusions/arriba/`, `results/fusions/fusioncatcher/`

---

### 4.9 CNV Inference

| Process | Tool | Description |
|---|---|---|
| `INFERCNV` | R / [InferCNV](https://github.com/broadinstitute/inferCNV) | Infer copy number variations from gene expression |

**Outputs:** `results/cnv_inference/`

---

### 4.10 WGCNA Co-expression Networks

| Process | Tool | Description |
|---|---|---|
| `WGCNA_ANALYSIS` | R / [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) | Weighted gene co-expression network analysis; module-trait correlations with ancestry and clinical variables |

**Parameters:** `wgcna_soft_power` (auto-detect if null), `wgcna_min_module_size` (30), `wgcna_merge_cut_height` (0.25), `wgcna_network_type` (signed)

**Outputs:** `results/wgcna/`

---

### 4.11 Pathway Enrichment

| Process | Tool | Description |
|---|---|---|
| `PATHWAY_ENRICHMENT` | R / [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/), [fgsea](https://bioconductor.org/packages/fgsea/) | Over-representation analysis (ORA) and gene set enrichment analysis (GSEA) |
| `GSVA_ANALYSIS` | R / [GSVA](https://bioconductor.org/packages/GSVA/) | Gene Set Variation Analysis -- per-sample pathway scores |

**Databases:** GO (BP, MF), KEGG, Reactome, MSigDB Hallmark, ImmuneSigDB

**Outputs:** `results/pathway_enrichment/ora_gsea/`, `results/pathway_enrichment/gsva/`

---

### 4.12 Immune Analysis

| Process | Tool | Description |
|---|---|---|
| `IMMUNE_DECONVOLUTION` | R / [immunedeconv](https://github.com/omnideconv/immunedeconv) | Multi-method immune cell deconvolution |
| `ESTIMATE_SCORES` | R / [ESTIMATE](https://bioinformatics.mdanderson.org/estimate/) | Tumor purity and immune/stromal scores from expression |

**Deconvolution methods:** [xCell](https://xcell.ucsf.edu/), [MCP-counter](https://github.com/ebecht/MCPcounter), [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/), [CIBERSORTx](https://cibersortx.stanford.edu/), [TIMER](http://timer.cistrome.org/)

All immune scores are corrected for tumor purity when `correct_tumor_purity = true` (default). ESTIMATE provides independent purity estimates from expression data.

**Outputs:** `results/immune_analysis/deconvolution/`, `results/immune_analysis/estimate/`

---

### 4.13 HLA Typing

| Process | Tool | Description |
|---|---|---|
| `ARCASHLA_EXTRACT` | [arcasHLA](https://github.com/RabadanLab/arcasHLA) 0.6.0 | Extract HLA reads from BAM |
| `ARCASHLA_GENOTYPE` | [arcasHLA](https://github.com/RabadanLab/arcasHLA) 0.6.0 | HLA class I and II genotyping (A, B, C, DPB1, DQB1, DQA1, DRB1) |
| `ARCASHLA_MERGE` | [arcasHLA](https://github.com/RabadanLab/arcasHLA) 0.6.0 | Merge per-sample HLA calls |
| `OPTITYPE` | [OptiType](https://github.com/FRED-2/OptiType) 1.3.5 | HLA class I genotyping (high precision) |

**Outputs:** `results/hla_typing/arcashla/`, `results/hla_typing/optitype/`

---

### 4.14 Neoantigen Prediction

All neoantigen prediction is performed from RNA-seq data alone -- no matched WGS/WES is required.

| Process | Tool | Description |
|---|---|---|
| `PVACSEQ` | [pVACseq](https://pvactools.readthedocs.io/) 4.2.0 | SNV/indel-derived neoantigen prediction with MHC binding |
| `NEOFUSE` | [NeoFuse](https://github.com/icbi-lab/NeoFuse) 1.0 | Fusion-derived neoantigen prediction |
| `SNAF_SPLICING_NEOANTIGENS` | [SNAF](https://github.com/frankligy/SNAF) | Splicing-derived neoantigen prediction |
| `MERGE_NEOANTIGENS` | Python (custom) | Merge neoantigens from all sources |

**Binding prediction algorithms:** MHCflurry, NetMHCpan, NetMHCIIpan

**Outputs:** `results/neoantigen/pvactools/`, `results/neoantigen/neofuse/`, `results/neoantigen/snaf/`, `results/neoantigen/merged/`

---

### 4.15 TCR/BCR Repertoire

| Process | Tool | Description |
|---|---|---|
| `TRUST4` | [TRUST4](https://github.com/liulab-dfci/TRUST4) 1.1.0 | TCR and BCR repertoire reconstruction from RNA-seq |
| `MERGE_TCR_REPORTS` | Python (custom) | Aggregate repertoire reports across samples |

**Outputs:** `results/tcr_repertoire/`

---

### 4.16 Molecular Subtyping

| Process | Tool | Description |
|---|---|---|
| `MOLECULAR_SUBTYPING` | R (custom) | Gene signature-based classification of disease subtypes |

For B-ALL (default), the module classifies samples into 12+ cytomolecular subtypes using curated gene signatures from the literature:

ETV6-RUNX1, BCR-ABL1, KMT2A rearranged, Hyperdiploid, Hypodiploid, Ph-like, DUX4, iAMP21, TCF3-PBX1, PAX5 alterations, MEF2D, ZNF384, NUTM1, and others.

**Outputs:** `results/molecular_subtyping/`

---

### 4.17 Sensitivity Analysis

| Process | Tool | Description |
|---|---|---|
| `SENSITIVITY_ANALYSIS` | R (custom) | Bootstrap resampling, timepoint stratification, ancestry threshold sensitivity |

Tests the robustness of DE results across:
- Bootstrap iterations (default: 1000)
- Timepoint stratification (diagnostic vs relapse)
- Varying ancestry proportion thresholds

**Outputs:** `results/sensitivity_analysis/`

---

### 4.18 Pharmacogenomics

| Process | Tool | Description |
|---|---|---|
| `PHARMACOGENOMICS` | R / [oncoPredict](https://github.com/Jing-Tao/oncoPredict), [DGIdb](https://www.dgidb.org/) | Drug sensitivity prediction from expression; drug-gene interaction mining |

**Drug response databases:** GDSC, CCLE, PRISM

**Outputs:** `results/pharmacogenomics/`

---

### 4.19 Validation

| Process | Tool | Description |
|---|---|---|
| `VALIDATE_DE_CONCORDANCE` | R (custom) | DESeq2 vs limma-voom concordance (overlap, rank correlation) |
| `VALIDATE_SEX_CHECK` | R (custom) | Verify reported sex against XIST/Y-chromosome gene expression |
| `VALIDATE_EXPRESSION_OUTLIERS` | R (custom) | PCA-based outlier detection (configurable SD threshold) |
| `VALIDATE_GENOMIC_INFLATION` | R (custom) | Genomic inflation factor (lambda) and QQ plots |
| `COMPILE_VALIDATION_REPORT` | Python (custom) | Aggregate validation results into a single report |

**Outputs:** `results/validation/de_concordance/`, `results/validation/sex_check/`, `results/validation/expression_outliers/`, `results/validation/genomic_inflation/`

---

### 4.20 Visualization and Reporting

| Process | Tool | Description |
|---|---|---|
| `GENERATE_REPORT` | Python (custom) | Integrated findings report across all modules |
| `VISUALIZATION_SUMMARY` | R (custom) | Integrated figures with figure index |
| `COMPILE_PIPELINE_LOG` | Python (custom) | Compile pipeline findings.md log |

**Outputs:** `results/report/`, `results/figures/`, `results/pipeline_log/`, `results/pipeline_info/`

---

## 5. Samplesheet Format

The input samplesheet is a CSV file with the following columns. Only `sample_id`, `bam`, and `bai` are required; all other columns are optional covariates.

| Column | Required | Type | Description |
|---|---|---|---|
| `sample_id` | Yes | string | Unique sample identifier |
| `bam` | Yes | path | Path to STAR-aligned, sorted, duplicate-marked BAM |
| `bai` | Yes | path | Path to corresponding BAM index |
| `batch` | No | string | Sequencing batch (for batch correction and DE covariate) |
| `sex` | No | string | Biological sex (M/F) |
| `age` | No | float | Age at sample collection |
| `tumor_purity` | No | float | Tumor purity estimate (0-1), e.g., from pathology or prior analysis |
| `blast_percentage` | No | float | Blast percentage from flow cytometry (B-ALL specific) |
| `cytomolecular_subgroup` | No | string | Known cytomolecular subtype (e.g., ETV6-RUNX1, Hyperdiploid) |
| `relapse_status` | No | string | Relapse status (relapse / no_relapse) |
| `adi_quartile` | No | string | Area Deprivation Index quartile (Q1-Q4) |
| `timepoint` | No | string | Sample timepoint (diagnostic / relapse) |
| `disease_stage` | No | string | Risk stratification (standard_risk / high_risk / very_high_risk) |

### Example

```csv
sample_id,bam,bai,batch,sex,age,tumor_purity,blast_percentage,cytomolecular_subgroup,relapse_status,adi_quartile,timepoint,disease_stage
SJBALL001,/path/to/SJBALL001.sorted.markdup.bam,/path/to/SJBALL001.sorted.markdup.bai,batch1,M,5,0.85,92,ETV6-RUNX1,no_relapse,Q2,diagnostic,standard_risk
SJBALL002,/path/to/SJBALL002.sorted.markdup.bam,/path/to/SJBALL002.sorted.markdup.bai,batch1,F,8,0.72,87,Hyperdiploid,relapse,Q1,diagnostic,standard_risk
SJBALL003,/path/to/SJBALL003.sorted.markdup.bam,/path/to/SJBALL003.sorted.markdup.bai,batch2,M,12,0.91,95,BCR-ABL1,no_relapse,Q3,diagnostic,high_risk
SJBALL004,/path/to/SJBALL004.sorted.markdup.bam,/path/to/SJBALL004.sorted.markdup.bai,batch2,F,3,0.78,76,KMT2A,relapse,Q1,relapse,high_risk
SJBALL005,/path/to/SJBALL005.sorted.markdup.bam,/path/to/SJBALL005.sorted.markdup.bai,batch1,M,7,0.65,88,Hypodiploid,no_relapse,Q4,diagnostic,very_high_risk
```

---

## 6. Contrast Definitions

Differential expression contrasts are defined via a JSON file passed with `--de_contrasts`. If no file is provided, the pipeline uses `assets/default_contrasts.json`.

### Contrast Types

| Type | Description | Fields |
|---|---|---|
| `categorical` | Two-group comparison | `variable`, `reference`, `target` |
| `continuous` | Regression on a continuous variable | `variable` |
| `categorical_multi` | One-vs-reference for each level of a multi-level factor | `variable`, `reference` |

### Example JSON

```json
[
    {
        "name": "relapse_vs_no_relapse",
        "variable": "relapse_status",
        "type": "categorical",
        "reference": "no_relapse",
        "target": "relapse",
        "description": "Differential expression: relapse vs. non-relapse samples"
    },
    {
        "name": "adi_q1_vs_rest",
        "variable": "adi_quartile",
        "type": "categorical",
        "reference": "rest",
        "target": "Q1",
        "description": "Area Deprivation Index lowest quartile vs. all others"
    },
    {
        "name": "pct_african_continuous",
        "variable": "pct_african",
        "type": "continuous",
        "description": "Continuous African ancestry proportion (from GRAF/ancestry inference)"
    },
    {
        "name": "graf_category_vs_EUR",
        "variable": "graf_category",
        "type": "categorical_multi",
        "reference": "EUR",
        "description": "GRAF ancestry categories (AFR_AM, LA1, LA2, EAS, SAS) each vs. EUR reference"
    },
    {
        "name": "timepoint_relapse_vs_diagnostic",
        "variable": "timepoint",
        "type": "categorical",
        "reference": "diagnostic",
        "target": "relapse",
        "description": "Relapse timepoint vs. diagnostic samples (paired when available)"
    }
]
```

### How contrasts work

- **`categorical`**: Fits a standard two-group comparison (target vs reference) controlling for covariates.
- **`continuous`**: Fits the variable as a continuous predictor in the DE model (e.g., ancestry proportions). Reports genes whose expression scales with the variable.
- **`categorical_multi`**: Automatically generates one contrast per non-reference level. For example, `graf_category_vs_EUR` with reference `EUR` creates separate contrasts for AFR_AM vs EUR, LA1 vs EUR, LA2 vs EUR, EAS vs EUR, and SAS vs EUR.

---

## 7. Covariates

The differential expression model includes configurable covariates specified via `--de_covariates` (comma-separated).

### Default DE model covariates

```
--de_covariates 'batch,sex,age,blast_percentage,tumor_purity,timepoint'
```

| Covariate | Type | Description |
|---|---|---|
| `batch` | Categorical | Sequencing batch -- accounts for technical variation |
| `sex` | Categorical | Biological sex (M/F) |
| `age` | Continuous | Age at sample collection |
| `blast_percentage` | Continuous | Leukemic blast percentage from clinical flow cytometry (B-ALL specific) |
| `tumor_purity` | Continuous | Computational tumor purity estimate (e.g., from ESTIMATE or pathology) |
| `timepoint` | Categorical | Sample collection timepoint (diagnostic/relapse) |

### blast_percentage vs tumor_purity

These two covariates measure different aspects and can both be included in the model:

- **`blast_percentage`** is a clinical measurement from flow cytometry that quantifies the proportion of leukemic blasts in the bone marrow aspirate. This is specific to hematologic malignancies.
- **`tumor_purity`** is a computational estimate of overall tumor content derived from expression data (via ESTIMATE) or pathology review. It captures broader tumor cellularity.

For non-hematologic cancers, drop `blast_percentage` and keep `tumor_purity`.

### Ancestry covariates

Ancestry is not specified in `--de_covariates` directly. Instead, ancestry proportions and categories from the GENOMICS subworkflow are automatically joined to the metadata and included in the DE model when ancestry-related contrasts are defined in the contrasts JSON.

---

## 8. Differential Variables

The pipeline supports the following categories of differential testing variables, each defined as a contrast in the JSON file.

### 8.1 Continuous Ancestry Proportions

Five continuous ancestry proportions derived from GRAF-pop SNP analysis:

| Variable | Description |
|---|---|
| `pct_african` | African ancestry proportion |
| `pct_amerindigenous` | Amerindigenous ancestry proportion |
| `pct_east_asian` | East Asian ancestry proportion |
| `pct_south_asian` | South Asian ancestry proportion |
| `pct_european` | European ancestry proportion |

These are used as continuous predictors in the DE model to identify genes whose expression varies linearly with genetic ancestry.

### 8.2 Categorical GRAF Ancestry Categories

Six discrete ancestry categories from supervised GRAF classification:

**EUR** (European), **AFR_AM** (African American), **LA1** (Latin American 1 -- higher European admixture), **LA2** (Latin American 2 -- higher Amerindigenous admixture), **EAS** (East Asian), **SAS** (South Asian)

Using `categorical_multi` contrast type, each non-reference category is compared against a reference (default: EUR).

### 8.3 Relapse Status

- **Contrast:** relapse vs no_relapse
- **Type:** categorical
- Identifies genes associated with treatment resistance and disease recurrence.

### 8.4 Area Deprivation Index (ADI)

- **Contrast:** Q1 (least deprived) vs rest
- **Type:** categorical
- Tests for gene expression differences associated with neighborhood socioeconomic deprivation.

### 8.5 Timepoint

- **Contrast:** relapse vs diagnostic
- **Type:** categorical (paired when longitudinal samples available)
- Identifies expression changes between disease stages.

### 8.6 Within Cytomolecular Subgroup Analyses

When `--cytomolecular_subgroups` is specified, the pipeline repeats all contrasts within each subgroup separately. This tests whether ancestry or clinical variable effects are consistent across molecular subtypes or are subtype-specific.

```bash
--cytomolecular_subgroups 'ETV6-RUNX1,Hyperdiploid,BCR-ABL1,KMT2A,Ph-like'
```

---

## 9. Parameters

### Input/Output

| Parameter | Default | Description |
|---|---|---|
| `--input` | (required) | Path to samplesheet CSV |
| `--outdir` | `./results` | Output directory |

### Reference Genome

| Parameter | Default | Description |
|---|---|---|
| `--genome` | `GRCh38` | Genome name |
| `--fasta` | (required) | Reference FASTA |
| `--fasta_fai` | null | FASTA index |
| `--gtf` | (required) | Gene annotation GTF |
| `--gene_bed` | null | Gene regions BED |
| `--star_index` | null | STAR genome index |
| `--known_snps` | null | dbSNP VCF for variant calling |
| `--known_snps_tbi` | null | dbSNP VCF index |

### Genome Build Detection and Liftover

| Parameter | Default | Description |
|---|---|---|
| `--target_genome_build` | `hg38` | Target build for analysis |
| `--chain_file` | null | Chain file for CrossMap liftover (e.g., hg19ToHg38.over.chain.gz) |
| `--auto_detect_build` | `true` | Auto-detect BAM genome build from header |
| `--force_liftover` | `false` | Force liftover even if builds match |

### Module Toggles

All modules can be individually enabled or disabled:

| Parameter | Default | Description |
|---|---|---|
| `--run_qc` | `true` | Quality control |
| `--run_ancestry` | `true` | Genetic ancestry inference |
| `--run_counting` | `true` | Read counting and normalization |
| `--run_variant_calling` | `true` | GATK variant calling |
| `--run_de` | `true` | Differential expression |
| `--run_splicing` | `true` | Differential splicing |
| `--run_fusions` | `true` | Gene fusion detection |
| `--run_cnv_inference` | `true` | CNV inference |
| `--run_wgcna` | `true` | WGCNA co-expression |
| `--run_pathway` | `true` | Pathway enrichment |
| `--run_immune` | `true` | Immune deconvolution |
| `--run_hla` | `true` | HLA typing |
| `--run_neoantigen` | `true` | Neoantigen prediction |
| `--run_tcr` | `true` | TCR/BCR repertoire |
| `--run_sensitivity` | `true` | Sensitivity analysis |
| `--run_pharmacogenomics` | `true` | Pharmacogenomics |
| `--run_molecular_subtyping` | `true` | Molecular subtyping |
| `--run_validation` | `true` | Validation checks |
| `--run_reporting` | `true` | Report generation |

### QC

| Parameter | Default | Description |
|---|---|---|
| `--skip_rseqc` | `false` | Skip RSeQC modules |
| `--skip_picard` | `false` | Skip Picard RNA-seq metrics |

### Counting

| Parameter | Default | Description |
|---|---|---|
| `--fc_extra_attributes` | `gene_name` | Extra GTF attributes for featureCounts |
| `--fc_group_features` | `gene_id` | GTF attribute for grouping |
| `--fc_count_type` | `exon` | Feature type to count |
| `--fc_strandedness` | `2` | 0=unstranded, 1=stranded, 2=reversely stranded |
| `--min_mapped_reads` | `500000` | Minimum mapped reads per sample |
| `--min_gene_counts` | `10` | Minimum counts across samples to keep a gene |
| `--min_samples_expressing` | `3` | Minimum samples with counts > 0 |
| `--normalization_method` | `vst` | Normalization: vst, rlog, tmm, tpm |

### Differential Expression

| Parameter | Default | Description |
|---|---|---|
| `--de_tool` | `both` | DE engine: deseq2, limma, both |
| `--padj_threshold` | `0.05` | Adjusted p-value cutoff |
| `--lfc_threshold` | `0.585` | log2 fold-change cutoff (log2(1.5)) |
| `--de_covariates` | `batch,sex,age,blast_percentage,tumor_purity,timepoint` | Model covariates |
| `--de_contrasts` | null | Path to contrasts JSON |
| `--cytomolecular_subgroups` | null | Subgroups for within-group analyses |

### Ancestry

| Parameter | Default | Description |
|---|---|---|
| `--graf_snp_bed` | null | GRAF ancestry-informative SNP BED positions |
| `--grafanc_data` | null | GRAF-anc data directory |
| `--ancestry_reference_panel` | null | Reference panel VCF |
| `--ancestry_reference_labels` | null | Population labels |
| `--graf_ancestry_categories` | `EUR,AFR_AM,LA1,LA2,EAS,SAS` | Categories to assign |
| `--ancestry_continuous_vars` | `pct_african,pct_amerindigenous,pct_south_asian,pct_east_asian,pct_european` | Continuous proportion variables |
| `--somalier_sites` | null | Somalier sites VCF |

### Splicing

| Parameter | Default | Description |
|---|---|---|
| `--splicing_tool` | `all` | rmats, leafcutter, spladder, bisbee, both, all |
| `--rmats_read_length` | `150` | Read length for rMATS |
| `--rmats_novel_ss` | `true` | Allow novel splice sites |
| `--leafcutter_min_coverage` | `20` | Minimum junction coverage |
| `--spladder_confidence` | `3` | SplAdder confidence (1-3) |
| `--bisbee_outlier_fdr` | `0.05` | FDR threshold for Bisbee outliers |

### Fusions

| Parameter | Default | Description |
|---|---|---|
| `--fusion_tool` | `both` | fusioncatcher, arriba, both |
| `--fusioncatcher_data` | null | FusionCatcher database directory |
| `--arriba_blacklist` | null | Arriba blacklist file |
| `--arriba_known_fusions` | null | Known fusions for Arriba |
| `--min_fusion_reads` | `3` | Minimum supporting reads |

### Variant Calling

| Parameter | Default | Description |
|---|---|---|
| `--variant_caller` | `gatk` | Variant caller (GATK HaplotypeCaller) |
| `--min_base_quality` | `20` | Minimum base quality |
| `--min_mapping_quality` | `20` | Minimum mapping quality |
| `--dbsnp` | null | dbSNP VCF for annotation |
| `--intervals` | null | Target regions BED |

### WGCNA

| Parameter | Default | Description |
|---|---|---|
| `--wgcna_min_module_size` | `30` | Minimum module size |
| `--wgcna_merge_cut_height` | `0.25` | Module merge cut height |
| `--wgcna_soft_power` | null | Soft-thresholding power (auto-detect if null) |
| `--wgcna_network_type` | `signed` | Network type |

### Pathway Enrichment

| Parameter | Default | Description |
|---|---|---|
| `--pathway_databases` | `GO_BP,GO_MF,KEGG,REACTOME,HALLMARK,IMMUNESIGDB` | MSigDB collections |
| `--msigdb_species` | `Homo sapiens` | Species for MSigDB |
| `--gsea_min_size` | `15` | Minimum gene set size for GSEA |
| `--gsea_max_size` | `500` | Maximum gene set size for GSEA |
| `--run_gsva` | `true` | Run GSVA analysis |

### Immune Analysis

| Parameter | Default | Description |
|---|---|---|
| `--immune_deconv_methods` | `cibersortx,xcell,mcpcounter,estimate,epic,timer` | Deconvolution methods |
| `--cibersortx_token` | null | CIBERSORTx API token |
| `--cibersortx_sigmatrix` | `LM22` | Signature matrix |
| `--correct_tumor_purity` | `true` | Adjust immune scores for tumor purity |

### HLA Typing

| Parameter | Default | Description |
|---|---|---|
| `--hla_tool` | `arcashla` | arcashla, optitype, both |
| `--arcashla_genes` | `A,B,C,DPB1,DQB1,DQA1,DRB1` | HLA genes to type |

### Neoantigen Prediction

| Parameter | Default | Description |
|---|---|---|
| `--neoantigen_tool` | `all` | pvacseq, neofuse, snaf, both, all |
| `--pvac_algorithms` | `MHCflurry,NetMHCpan,NetMHCIIpan` | MHC binding predictors |
| `--binding_threshold` | `500` | Binding affinity cutoff (nM) |
| `--percentile_threshold` | `2.0` | Binding percentile cutoff |
| `--min_variant_depth` | `10` | Minimum variant read depth |
| `--min_vaf` | `0.05` | Minimum variant allele frequency |
| `--snaf_db` | null | SNAF reference database directory |

### Pharmacogenomics

| Parameter | Default | Description |
|---|---|---|
| `--drug_response_db` | `GDSC` | Drug response database: GDSC, CCLE, PRISM, all |
| `--dgidb_interactions` | `true` | Query DGIdb for drug-gene interactions |

### Sensitivity Analysis

| Parameter | Default | Description |
|---|---|---|
| `--sensitivity_variables` | `timepoint,ancestry_proportion` | Variables to test |
| `--bootstrap_iterations` | `1000` | Number of bootstrap resamples |

### Batch Correction

| Parameter | Default | Description |
|---|---|---|
| `--batch_correction` | `combat_seq` | Method: combat_seq, sva, none |
| `--batch_variable` | `batch` | Column name for batch variable |

### Validation

| Parameter | Default | Description |
|---|---|---|
| `--outlier_sd_threshold` | `3` | SD threshold for PCA outlier detection |
| `--min_sample_correlation` | `0.8` | Minimum pairwise sample correlation |

### Resource Limits

| Parameter | Default | Description |
|---|---|---|
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_time` | `240.h` | Maximum wall time per process |

---

## 10. Output Directory Structure

```
results/
|-- qc/
|   |-- samtools/                  # flagstat, idxstats, stats per sample
|   |-- rseqc/                     # bam_stat, read_distribution, infer_experiment
|   |-- picard/                    # CollectRnaSeqMetrics per sample
|   `-- multiqc/                   # Aggregated MultiQC HTML report
|
|-- counting/
|   |-- featurecounts/             # Per-sample gene counts + summary
|   |-- normalized/                # VST/rlog/TMM/TPM normalized count matrix
|   `-- batch_corrected/           # ComBat-seq corrected counts
|
|-- ancestry/
|   |-- graf_snps/                 # Per-sample allele counts at GRAF SNPs
|   |-- grafanc/                   # GRAF-anc raw output
|   |-- proportions/               # Ancestry proportions + categories TSV
|   `-- somalier/                  # Somalier relatedness + ancestry QC
|
|-- variant_calling/
|   |-- raw/                       # Per-sample raw VCFs (HaplotypeCaller)
|   |-- filtered/                  # Hard-filtered VCFs
|   `-- tmb/                       # TMB scores per sample
|
|-- differential_expression/
|   |-- deseq2/                    # DE results per contrast (tables + plots)
|   `-- limma/                     # Limma-voom results per contrast
|
|-- splicing/
|   |-- rmats/                     # Differential splicing events (SE, MXE, etc.)
|   |-- leafcutter/                # Intron cluster counts + differential usage
|   |-- spladder/
|   |   |-- graphs/                # Augmented splicing graphs
|   |   `-- testing/               # Statistical testing results
|   `-- bisbee/
|       |-- differential/          # Beta-binomial differential splicing
|       |-- protein_effects/       # Predicted protein-level effects
|       `-- outliers/              # Per-sample outlier splicing events
|
|-- fusions/
|   |-- arriba/                    # Fusion calls + visualization PDFs
|   `-- fusioncatcher/             # Fusion calls + supporting evidence
|
|-- cnv_inference/                 # InferCNV heatmaps + CNV calls
|
|-- wgcna/                         # Module assignments, eigengenes, trait correlations
|
|-- pathway_enrichment/
|   |-- ora_gsea/                  # ORA + GSEA results per contrast per database
|   `-- gsva/                      # Per-sample GSVA pathway activity scores
|
|-- immune_analysis/
|   |-- deconvolution/             # Immune cell scores (xCell, MCP, EPIC, CIBERSORT, TIMER)
|   `-- estimate/                  # ESTIMATE purity, immune, stromal scores
|
|-- hla_typing/
|   |-- arcashla/                  # Per-sample + merged HLA genotypes
|   `-- optitype/                  # OptiType HLA class I calls
|
|-- neoantigen/
|   |-- pvactools/                 # SNV/indel neoantigens per sample
|   |-- neofuse/                   # Fusion-derived neoantigens
|   |-- snaf/                      # Splicing-derived neoantigens
|   `-- merged/                    # Combined neoantigen report
|
|-- tcr_repertoire/                # TRUST4 TCR/BCR reports + merged summary
|
|-- molecular_subtyping/           # Subtype classifications + confidence scores
|
|-- sensitivity_analysis/          # Bootstrap, timepoint, threshold sensitivity results
|
|-- pharmacogenomics/              # Drug sensitivity predictions + DGIdb interactions
|
|-- validation/
|   |-- de_concordance/            # DESeq2 vs limma overlap and rank correlation
|   |-- sex_check/                 # XIST/Y-gene expression vs reported sex
|   |-- expression_outliers/       # PCA outlier detection results
|   `-- genomic_inflation/         # Lambda values + QQ plots per contrast
|
|-- report/                        # Integrated findings report
|-- figures/                       # Publication-ready figures + figure index HTML
|-- pipeline_log/                  # Compiled findings.md log
|
`-- pipeline_info/
    |-- timeline_<timestamp>.html  # Nextflow execution timeline
    |-- report_<timestamp>.html    # Nextflow execution report
    |-- trace_<timestamp>.txt      # Process-level resource trace
    `-- dag_<timestamp>.html       # Pipeline DAG visualization
```

---

## 11. Configurability for Other Cancers

While the default configuration targets pediatric B-ALL, the pipeline is designed to be adapted for any cancer type. The main adjustments required are:

### 1. Covariates

Remove B-ALL specific covariates and add cancer-appropriate ones:

```bash
# Solid tumor example (no blast_percentage)
--de_covariates 'batch,sex,age,tumor_purity,tumor_stage'

# AML example
--de_covariates 'batch,sex,age,blast_percentage,tumor_purity,FAB_subtype'
```

### 2. Contrasts

Write a custom contrasts JSON for your disease-specific comparisons:

```json
[
    {
        "name": "tumor_vs_normal",
        "variable": "sample_type",
        "type": "categorical",
        "reference": "normal",
        "target": "tumor"
    },
    {
        "name": "stage_IV_vs_I",
        "variable": "tumor_stage",
        "type": "categorical",
        "reference": "I",
        "target": "IV"
    }
]
```

### 3. Molecular Subtyping

For non-B-ALL cancers, either:
- Disable the module: `--run_molecular_subtyping false`
- Provide custom gene signature files appropriate for your cancer type

### 4. Samplesheet Columns

Add any custom clinical columns to the samplesheet. They will be available as metadata fields. The pipeline does not enforce a fixed set of column names beyond `sample_id`, `bam`, and `bai`.

### 5. Immune Deconvolution

The immune deconvolution methods are cancer-agnostic. Adjust the ESTIMATE platform if not using Illumina:

```bash
--estimate_platform 'affymetrix'  # or 'illumina' (default)
```

### 6. Pathway Databases

Customize pathway databases for your disease context:

```bash
--pathway_databases 'GO_BP,KEGG,REACTOME,HALLMARK'
```

---

## 12. Profiles

The pipeline ships with the following execution profiles, selectable via `-profile`:

| Profile | Description | Usage |
|---|---|---|
| `docker` | Run all processes in Docker containers | `-profile docker` |
| `singularity` | Run all processes in Singularity containers (HPC compatible) | `-profile singularity` |
| `conda` | Use Conda environments | `-profile conda` |
| `slurm` | Submit jobs to a SLURM scheduler | `-profile singularity,slurm` |
| `lsf` | Submit jobs to an LSF scheduler | `-profile singularity,lsf` |
| `local` | Run on local machine | `-profile docker,local` |
| `test` | Minimal test data, reduced resources (CI/CD) | `-profile test,docker` |

### Combining profiles

Profiles are composable. Combine a container profile with a scheduler profile:

```bash
# Singularity on a SLURM cluster
nextflow run main.nf -profile singularity,slurm --input samplesheet.csv ...

# Docker locally
nextflow run main.nf -profile docker,local --input samplesheet.csv ...
```

### SLURM configuration

The SLURM profile (`conf/profiles/slurm.config`) routes high-memory jobs (>=64 GB) to the `bigmem` partition and standard jobs to the `normal` partition. Adjust `clusterOptions` for your site:

```groovy
process {
    executor = 'slurm'
    queue    = { task.memory >= 64.GB ? 'bigmem' : 'normal' }
    clusterOptions = '--account=your_account'
}
```

### Test profile

The test profile disables resource-intensive modules (WGCNA, CNV inference, neoantigen prediction, pharmacogenomics, sensitivity analysis, molecular subtyping) and reduces resource caps to 2 CPUs, 6 GB memory, and 6 hours max wall time.

### Resource labels

Processes use resource labels defined in `conf/base.config`:

| Label | CPUs | Memory | Time |
|---|---|---|---|
| `process_low` | 2 | 4 GB | 2 h |
| `process_medium` | 8 | 32 GB | 8 h |
| `process_high` | 16 | 64 GB | 24 h |
| `process_very_high` | 16 | 128 GB | 48 h |
| `process_long` | default | default | 96 h |
| `process_single` | 1 | 8 GB | default |
| `process_gpu` | 4 | 32 GB | 8 h |

All resources auto-retry up to 2 times on transient failures (exit codes 143, 137, 104, 134, 139, 140) with increasing resource allocation per attempt.

---

## 13. Pipeline Logging

### findings.md Tracking System

Each module that generates notable findings appends structured entries to a `findings.md` log file. This provides a human-readable audit trail of pipeline decisions and significant results.

The log captures:
- **QC flags**: Samples failing quality thresholds, low read counts, strandedness mismatches
- **Ancestry assignments**: Per-sample ancestry classification and SNP coverage statistics
- **Genome build**: Detected build per BAM and whether liftover was performed
- **DE summary**: Number of significant genes per contrast, concordance between DESeq2 and limma-voom
- **Splicing outliers**: Bisbee outlier events flagged per sample
- **Fusion calls**: High-confidence fusions detected by both callers
- **Neoantigen counts**: Number of predicted neoantigens per sample per source (SNV, fusion, splicing)
- **Validation warnings**: Sex mismatches, PCA outliers, high genomic inflation (lambda > 1.1)

The `COMPILE_PIPELINE_LOG` process aggregates all findings into a final report at `results/pipeline_log/`.

### Nextflow execution reports

The pipeline also generates standard Nextflow execution reports:

| Report | Location | Description |
|---|---|---|
| Timeline | `results/pipeline_info/timeline_*.html` | Visual timeline of all processes |
| Execution report | `results/pipeline_info/report_*.html` | Resource usage summary |
| Trace | `results/pipeline_info/trace_*.txt` | Detailed per-process resource trace |
| DAG | `results/pipeline_info/dag_*.html` | Directed acyclic graph of the pipeline |

---

## 14. Visualizations Produced

The pipeline generates the following figures and plots across its modules.

### QC and Preprocessing

- MultiQC aggregate report (HTML)
- Per-sample alignment statistics bar plots
- Read distribution pie charts (CDS, UTR, intron, intergenic)
- Strandedness inference plots

### Ancestry

- PCA scatter plots of ancestry-informative SNPs (colored by GRAF category)
- Ancestry proportion stacked bar charts per sample
- Somalier relatedness heatmap

### Differential Expression

- Volcano plots per contrast (DESeq2 and limma-voom)
- MA plots per contrast
- Heatmaps of top DE genes
- PCA plots colored by contrast variable
- Venn diagram of DESeq2 vs limma-voom overlap
- Upset plots of gene overlap across contrasts

### Splicing

- Sashimi plots for significant differential splicing events
- Leafcutter cluster plots
- Bisbee outlier score distributions

### Fusions

- Arriba fusion circos plots per sample
- Fusion waterfall plots (recurrence across samples)

### CNV Inference

- InferCNV heatmaps (genome-wide expression-derived CNV)

### WGCNA

- Scale-free topology fit plots (soft power selection)
- Module dendrogram with color assignments
- Module-trait correlation heatmap (including ancestry traits)
- Module eigengene bar plots

### Pathway Enrichment

- Dot plots for top enriched GO/KEGG/Reactome terms
- GSEA running enrichment score plots for top gene sets
- Hallmark and ImmuneSigDB enrichment bar plots
- GSVA heatmap of per-sample pathway activity scores

### Immune Analysis

- Immune cell composition stacked bar charts per sample
- Immune cell score heatmap across methods
- ESTIMATE score box plots (immune, stromal, purity)
- Immune score vs ancestry proportion scatter plots

### Neoantigen Prediction

- Neoantigen burden bar charts per sample
- Binding affinity distribution plots
- Neoantigen source breakdown (SNV vs fusion vs splicing)

### TCR/BCR Repertoire

- Clonotype diversity plots (Shannon, Simpson indices)
- V-gene usage bar charts
- CDR3 length distribution histograms

### Pharmacogenomics

- Predicted drug sensitivity heatmap
- Drug-gene interaction network plots

### Validation

- DESeq2 vs limma-voom rank correlation scatter
- QQ plots with genomic inflation (lambda) per contrast
- Sex check scatter (XIST vs Y-chromosome expression)
- PCA outlier detection plots with SD threshold boundaries

### Reporting

- Figure index HTML page with thumbnails and links to all generated figures

---

## 15. Credits and References

### Pipeline

This pipeline was developed by the Cancer Genomics Lab.

### Core Tools

| Tool | Publication / Reference |
|---|---|
| [Nextflow](https://www.nextflow.io/) | Di Tommaso et al. (2017) *Nature Biotechnology* 35:316-319 |
| [samtools](http://www.htslib.org/) | Danecek et al. (2021) *GigaScience* 10:giab008 |
| [RSeQC](https://rseqc.sourceforge.net/) | Wang et al. (2012) *Bioinformatics* 28:2184-2185 |
| [Picard](https://broadinstitute.github.io/picard/) | Broad Institute |
| [MultiQC](https://multiqc.info/) | Ewels et al. (2016) *Bioinformatics* 32:3047-3048 |
| [Subread/featureCounts](https://subread.sourceforge.net/) | Liao et al. (2014) *Bioinformatics* 30:923-930 |
| [DESeq2](https://bioconductor.org/packages/DESeq2/) | Love et al. (2014) *Genome Biology* 15:550 |
| [limma](https://bioconductor.org/packages/limma/) | Ritchie et al. (2015) *Nucleic Acids Research* 43:e47 |
| [ComBat-seq](https://bioconductor.org/packages/sva/) | Zhang et al. (2020) *NAR Genomics and Bioinformatics* 2:lqaa078 |
| [GATK](https://gatk.broadinstitute.org/) | McKenna et al. (2010) *Genome Research* 20:1297-1303 |
| [CrossMap](https://crossmap.sourceforge.net/) | Zhao et al. (2014) *Bioinformatics* 30:1006-1007 |

### Ancestry and Population Genetics

| Tool | Publication / Reference |
|---|---|
| [GRAF-pop](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/Software.cgi) | Jin et al. (2019) *Genome Research* 29:1622-1638 |
| [Somalier](https://github.com/brentp/somalier) | Pedersen et al. (2020) *Genome Medicine* 12:62 |

### Splicing

| Tool | Publication / Reference |
|---|---|
| [rMATS](http://rnaseq-mats.sourceforge.net/) | Shen et al. (2014) *PNAS* 111:E5593-E5601 |
| [Leafcutter](https://davidaknowles.github.io/leafcutter/) | Li et al. (2018) *Nature Genetics* 50:151-158 |
| [SplAdder](https://github.com/ratschlab/spladder) | Kahles et al. (2016) *Bioinformatics* 32:1840-1847 |
| [Bisbee](https://github.com/tgen/bisbee) | Szabo et al. (2023) *Bioinformatics* |

### Fusions

| Tool | Publication / Reference |
|---|---|
| [Arriba](https://arriba.readthedocs.io/) | Uhrig et al. (2021) *Genome Research* 31:448-460 |
| [FusionCatcher](https://github.com/ndaniel/fusioncatcher) | Nicorici et al. (2014) *bioRxiv* |

### Immunogenomics

| Tool | Publication / Reference |
|---|---|
| [immunedeconv](https://github.com/omnideconv/immunedeconv) | Sturm et al. (2019) *Bioinformatics* 35:i436-i445 |
| [ESTIMATE](https://bioinformatics.mdanderson.org/estimate/) | Yoshihara et al. (2013) *Nature Communications* 4:2612 |
| [xCell](https://xcell.ucsf.edu/) | Aran et al. (2017) *Genome Biology* 18:220 |
| [CIBERSORTx](https://cibersortx.stanford.edu/) | Newman et al. (2019) *Nature Biotechnology* 37:711-718 |
| [MCP-counter](https://github.com/ebecht/MCPcounter) | Becht et al. (2016) *Genome Biology* 17:218 |
| [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) | Racle et al. (2017) *eLife* 6:e26476 |
| [TIMER](http://timer.cistrome.org/) | Li et al. (2017) *Genome Biology* 18:234 |

### HLA and Neoantigens

| Tool | Publication / Reference |
|---|---|
| [arcasHLA](https://github.com/RabadanLab/arcasHLA) | Orenbuch et al. (2020) *Bioinformatics* 36:33-40 |
| [OptiType](https://github.com/FRED-2/OptiType) | Szolek et al. (2014) *Bioinformatics* 30:3310-3316 |
| [pVACseq](https://pvactools.readthedocs.io/) | Hundal et al. (2020) *Cancer Immunology Research* 8:409-420 |
| [NeoFuse](https://github.com/icbi-lab/NeoFuse) | Fang et al. (2022) *Bioinformatics* 38:3983-3985 |
| [SNAF](https://github.com/frankligy/SNAF) | Li & Bhatt (2023) *Science Advances* 9:eade2886 |

### TCR/BCR and Co-expression

| Tool | Publication / Reference |
|---|---|
| [TRUST4](https://github.com/liulab-dfci/TRUST4) | Song et al. (2021) *Nature Methods* 18:820-823 |
| [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) | Langfelder & Horvath (2008) *BMC Bioinformatics* 9:559 |

### Pathway and Gene Set Analysis

| Tool | Publication / Reference |
|---|---|
| [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) | Wu et al. (2021) *The Innovation* 2:100141 |
| [fgsea](https://bioconductor.org/packages/fgsea/) | Korotkevich et al. (2021) *bioRxiv* |
| [GSVA](https://bioconductor.org/packages/GSVA/) | Hanzelmann et al. (2013) *BMC Bioinformatics* 14:7 |
| [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) | Liberzon et al. (2015) *Cell Systems* 1:417-425 |

### CNV and Pharmacogenomics

| Tool | Publication / Reference |
|---|---|
| [InferCNV](https://github.com/broadinstitute/inferCNV) | Patel et al. (2014) *Science* 344:1396-1401 |
| [oncoPredict](https://github.com/Jing-Tao/oncoPredict) | Maeser et al. (2021) *Briefings in Bioinformatics* 22:bbab260 |
| [DGIdb](https://www.dgidb.org/) | Freshour et al. (2021) *Nucleic Acids Research* 49:D1144-D1151 |

---

## License

Please refer to the LICENSE file in this repository.

## Issues and Contributions

Report issues and feature requests via the repository issue tracker. Contributions via pull requests are welcome.
