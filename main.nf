#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Cancer Bulk RNA-Seq Analysis Pipeline with Ancestry Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Modules:
      1.  QC & BAM metrics
      2.  Genetic ancestry inference (Somalier + GRAF-pop SNPs)
      3.  Read counting (Subread featureCounts) & normalization
      4.  Batch effect correction (ComBat-seq / SVA)
      5.  Variant calling from RNA-seq (GATK best practices)
      6.  Differential expression (DESeq2 + limma-voom)
      7.  Differential RNA splicing (rMATS / Leafcutter)
      8.  Gene fusion detection (FusionCatcher + Arriba)
      9.  CNV inference from RNA-seq (InferCNV)
      10. WGCNA co-expression networks
      11. Pathway enrichment (ORA, GSEA, GSVA)
      12. Immune deconvolution & microenvironment
      13. HLA typing (arcasHLA / OptiType)
      14. Neoantigen prediction (pVACseq / NeoFuse)
      15. TCR/BCR repertoire (TRUST4)
      16. Molecular subtyping / cell-of-origin
      17. Sensitivity analysis
      18. Pharmacogenomics & drug target prediction
      19. TMB estimation
      20. Integrated report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESSING        } from './subworkflows/local/preprocessing'
include { GENOMICS             } from './subworkflows/local/genomics'
include { EXPRESSION_ANALYSIS  } from './subworkflows/local/expression_analysis'
include { STRUCTURAL_VARIANTS  } from './subworkflows/local/structural_variants'
include { IMMUNOGENOMICS       } from './subworkflows/local/immunogenomics'
include { CLINICAL_ANALYSIS    } from './subworkflows/local/clinical_analysis'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC           } from './modules/local/qc/main'
include { GENERATE_REPORT   } from './modules/local/reporting/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate required params
if (!params.input) { exit 1, "ERROR: --input samplesheet not specified" }
if (!params.fasta) { exit 1, "ERROR: --fasta reference genome not specified" }
if (!params.gtf)   { exit 1, "ERROR: --gtf gene annotation not specified" }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARSE SAMPLESHEET
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parse_samplesheet(samplesheet_path) {
    Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [:]
            meta.id                     = row.sample_id
            meta.batch                  = row.batch ?: 'NA'
            meta.sex                    = row.sex ?: 'NA'
            meta.age                    = row.age ? row.age.toFloat() : null
            meta.tumor_purity           = row.tumor_purity ? row.tumor_purity.toFloat() : null
            meta.blast_percentage       = row.blast_percentage ? row.blast_percentage.toFloat() : null
            meta.cytomolecular_subgroup = row.cytomolecular_subgroup ?: 'NA'
            meta.relapse_status         = row.relapse_status ?: 'NA'
            meta.adi_quartile           = row.adi_quartile ?: 'NA'
            meta.timepoint              = row.timepoint ?: 'diagnostic'
            meta.disease_stage          = row.disease_stage ?: 'NA'

            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)

            return [ meta, bam, bai ]
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // ---- Reference files ----
    ch_fasta          = Channel.fromPath(params.fasta, checkIfExists: true).collect()
    ch_fasta_fai      = params.fasta_fai ? Channel.fromPath(params.fasta_fai, checkIfExists: true).collect() : Channel.empty()
    ch_gtf            = Channel.fromPath(params.gtf, checkIfExists: true).collect()
    ch_known_snps     = params.known_snps ? Channel.fromPath(params.known_snps, checkIfExists: true).collect() : Channel.empty()
    ch_known_snps_tbi = params.known_snps_tbi ? Channel.fromPath(params.known_snps_tbi, checkIfExists: true).collect() : Channel.empty()
    ch_graf_snp_bed   = params.graf_snp_bed ? Channel.fromPath(params.graf_snp_bed, checkIfExists: true).collect() : Channel.empty()

    // ---- Parse samplesheet ----
    ch_samples = parse_samplesheet(params.input)

    // Separate meta and BAM channels
    ch_bam_bai = ch_samples.map { meta, bam, bai -> [ meta, bam, bai ] }

    // Collect all metadata for downstream analyses
    ch_metadata = ch_samples
        .map { meta, bam, bai -> meta }
        .collect()

    // ---- SUBWORKFLOW 1: Preprocessing (QC + Counting + Normalization) ----
    PREPROCESSING(
        ch_bam_bai,
        ch_gtf,
        ch_fasta
    )
    ch_counts_raw        = PREPROCESSING.out.counts_raw
    ch_counts_normalized = PREPROCESSING.out.counts_normalized
    ch_count_matrix      = PREPROCESSING.out.count_matrix
    ch_qc_reports        = PREPROCESSING.out.qc_reports

    // ---- SUBWORKFLOW 2: Genomics (Ancestry + Variant Calling + HLA) ----
    GENOMICS(
        ch_bam_bai,
        ch_fasta,
        ch_fasta_fai,
        ch_known_snps,
        ch_known_snps_tbi,
        ch_graf_snp_bed
    )
    ch_ancestry          = GENOMICS.out.ancestry_proportions
    ch_variants          = GENOMICS.out.filtered_vcf
    ch_hla_types         = GENOMICS.out.hla_types
    ch_tmb               = GENOMICS.out.tmb_scores

    // ---- SUBWORKFLOW 3: Expression Analysis (DE + WGCNA + Pathway + GSVA) ----
    EXPRESSION_ANALYSIS(
        ch_count_matrix,
        ch_counts_normalized,
        ch_metadata,
        ch_ancestry
    )
    ch_de_results        = EXPRESSION_ANALYSIS.out.de_results
    ch_wgcna_modules     = EXPRESSION_ANALYSIS.out.wgcna_modules
    ch_pathway_results   = EXPRESSION_ANALYSIS.out.pathway_results
    ch_gsva_scores       = EXPRESSION_ANALYSIS.out.gsva_scores

    // ---- SUBWORKFLOW 4: Structural Variants (Splicing + Fusions + CNV) ----
    STRUCTURAL_VARIANTS(
        ch_bam_bai,
        ch_fasta,
        ch_gtf,
        ch_metadata
    )
    ch_splicing          = STRUCTURAL_VARIANTS.out.splicing_results
    ch_fusions           = STRUCTURAL_VARIANTS.out.fusion_results
    ch_cnv               = STRUCTURAL_VARIANTS.out.cnv_results

    // ---- SUBWORKFLOW 5: Immunogenomics (Immune + Neoantigen + TCR) ----
    IMMUNOGENOMICS(
        ch_bam_bai,
        ch_count_matrix,
        ch_counts_normalized,
        ch_variants,
        ch_hla_types,
        ch_fusions,
        ch_metadata
    )
    ch_immune            = IMMUNOGENOMICS.out.immune_scores
    ch_neoantigens       = IMMUNOGENOMICS.out.neoantigens
    ch_tcr               = IMMUNOGENOMICS.out.tcr_repertoire

    // ---- SUBWORKFLOW 6: Clinical Analysis (Sensitivity + Pharma + Subtyping) ----
    CLINICAL_ANALYSIS(
        ch_count_matrix,
        ch_counts_normalized,
        ch_de_results,
        ch_ancestry,
        ch_immune,
        ch_metadata
    )
    ch_sensitivity       = CLINICAL_ANALYSIS.out.sensitivity_results
    ch_pharma            = CLINICAL_ANALYSIS.out.pharma_results
    ch_subtypes          = CLINICAL_ANALYSIS.out.subtype_results

    // ---- MultiQC ----
    if (params.run_reporting) {
        ch_multiqc_files = Channel.empty()
            .mix(ch_qc_reports.collect().ifEmpty([]))
        MULTIQC(ch_multiqc_files.collect())
    }

    // ---- Integrated Report ----
    if (params.run_reporting) {
        ch_all_results = Channel.empty()
            .mix(
                ch_de_results.ifEmpty([]),
                ch_pathway_results.ifEmpty([]),
                ch_immune.ifEmpty([]),
                ch_tmb.ifEmpty([]),
                ch_neoantigens.ifEmpty([]),
                ch_pharma.ifEmpty([]),
                ch_sensitivity.ifEmpty([])
            )
        GENERATE_REPORT(ch_all_results.collect())
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ON COMPLETE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    log.info ""
    log.info "========================================="
    log.info "Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration          : ${workflow.duration}"
    log.info "Output directory  : ${params.outdir}"
    log.info "========================================="
    log.info ""
}

workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
}
