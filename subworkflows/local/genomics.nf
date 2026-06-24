/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOMICS Subworkflow
    Ancestry inference (GRAF-anc + Somalier) + Variant calling + HLA typing + TMB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXTRACT_GRAF_SNPS      } from '../../modules/local/ancestry/main'
include { MERGE_GRAF_VCFS        } from '../../modules/local/ancestry/main'
include { GRAFANC_RUN            } from '../../modules/local/ancestry/main'
include { ANCESTRY_INFERENCE     } from '../../modules/local/ancestry/main'
include { SOMALIER_EXTRACT       } from '../../modules/local/ancestry/main'
include { SOMALIER_RELATE        } from '../../modules/local/ancestry/main'
include { GATK_SPLITNCIGARREADS  } from '../../modules/local/variant_calling/main'
include { GATK_BASERECALIBRATOR  } from '../../modules/local/variant_calling/main'
include { GATK_APPLYBQSR         } from '../../modules/local/variant_calling/main'
include { GATK_HAPLOTYPECALLER   } from '../../modules/local/variant_calling/main'
include { GATK_VARIANTFILTRATION } from '../../modules/local/variant_calling/main'
include { TMB_ESTIMATION         } from '../../modules/local/variant_calling/main'
include { CHASMPLUS              } from '../../modules/local/variant_calling/main'
include { CREATE_MAF             } from '../../modules/local/variant_calling/main'
include { ARCASHLA_EXTRACT       } from '../../modules/local/hla_typing/main'
include { ARCASHLA_GENOTYPE      } from '../../modules/local/hla_typing/main'
include { ARCASHLA_MERGE         } from '../../modules/local/hla_typing/main'
include { OPTITYPE               } from '../../modules/local/hla_typing/main'

workflow GENOMICS {

    take:
    ch_bam_bai        // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta          // channel: path(fasta)
    ch_fasta_fai      // channel: path(fasta_fai)
    ch_known_snps     // channel: path(known_snps)
    ch_known_snps_tbi // channel: path(known_snps_tbi)
    ch_graf_snp_bed   // channel: path(graf_snp_bed)

    main:
    ch_versions = Channel.empty()
    ch_findings = Channel.empty()

    // ========================================
    // ANCESTRY INFERENCE (GRAF-anc + Somalier)
    // ========================================
    if (params.run_ancestry) {
        // --- GRAF-anc: Extract 282K ancestry-informative SNPs ---
        EXTRACT_GRAF_SNPS(ch_bam_bai, ch_graf_snp_bed.first(), ch_fasta.first())
        ch_versions = ch_versions.mix(EXTRACT_GRAF_SNPS.out.versions.first())

        // Merge per-sample VCFs into multi-sample VCF for GRAF-anc
        MERGE_GRAF_VCFS(
            EXTRACT_GRAF_SNPS.out.vcf.map{ meta, vcf -> vcf }.collect(),
            EXTRACT_GRAF_SNPS.out.tbi.map{ meta, tbi -> tbi }.collect()
        )
        ch_versions = ch_versions.mix(MERGE_GRAF_VCFS.out.versions)

        // Run GRAF-anc on merged VCF
        if (params.grafanc_data) {
            GRAFANC_RUN(
                MERGE_GRAF_VCFS.out.merged_vcf,
                MERGE_GRAF_VCFS.out.merged_tbi,
                Channel.fromPath(params.grafanc_data).first()
            )
            ch_versions = ch_versions.mix(GRAFANC_RUN.out.versions)

            // Full ancestry inference with GRAF-anc results
            ANCESTRY_INFERENCE(
                GRAFANC_RUN.out.results,
                EXTRACT_GRAF_SNPS.out.allele_counts.map{ meta, f -> f }.collect(),
                params.ancestry_reference_panel ? Channel.fromPath(params.ancestry_reference_panel).first() : file('NO_REF_PANEL'),
                params.ancestry_reference_labels ? Channel.fromPath(params.ancestry_reference_labels).first() : file('NO_REF_LABELS')
            )
        } else {
            // Fallback: ancestry inference from allele counts only (no GRAF-anc binary)
            ANCESTRY_INFERENCE(
                file('NO_GRAFANC_RESULTS'),
                EXTRACT_GRAF_SNPS.out.allele_counts.map{ meta, f -> f }.collect(),
                params.ancestry_reference_panel ? Channel.fromPath(params.ancestry_reference_panel).first() : file('NO_REF_PANEL'),
                params.ancestry_reference_labels ? Channel.fromPath(params.ancestry_reference_labels).first() : file('NO_REF_LABELS')
            )
        }
        ch_versions = ch_versions.mix(ANCESTRY_INFERENCE.out.versions)
        ch_findings = ch_findings.mix(ANCESTRY_INFERENCE.out.findings)

        // --- Somalier: Sample QC / relatedness checking ---
        if (params.somalier_sites) {
            SOMALIER_EXTRACT(
                ch_bam_bai,
                ch_fasta.first(),
                Channel.fromPath(params.somalier_sites).first()
            )
            ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions.first())

            SOMALIER_RELATE(
                SOMALIER_EXTRACT.out.extracted.map{ meta, f -> f }.collect()
            )
            ch_versions = ch_versions.mix(SOMALIER_RELATE.out.versions)
        }
    }

    // ========================================
    // VARIANT CALLING (GATK RNA-seq best practices)
    // ========================================
    if (params.run_variant_calling) {
        // Step 1: Split N CIGAR reads
        GATK_SPLITNCIGARREADS(ch_bam_bai, ch_fasta.first(), ch_fasta_fai.first())
        ch_versions = ch_versions.mix(GATK_SPLITNCIGARREADS.out.versions.first())

        // Step 2: Base recalibration (if known sites available)
        if (params.known_snps) {
            GATK_BASERECALIBRATOR(
                GATK_SPLITNCIGARREADS.out.bam,
                ch_fasta.first(),
                ch_fasta_fai.first(),
                ch_known_snps.first(),
                ch_known_snps_tbi.first()
            )
            ch_versions = ch_versions.mix(GATK_BASERECALIBRATOR.out.versions.first())

            ch_bam_recal = GATK_SPLITNCIGARREADS.out.bam
                .join(GATK_BASERECALIBRATOR.out.table)

            GATK_APPLYBQSR(ch_bam_recal, ch_fasta.first(), ch_fasta_fai.first())
            ch_versions = ch_versions.mix(GATK_APPLYBQSR.out.versions.first())

            ch_bam_for_calling = GATK_APPLYBQSR.out.bam
        } else {
            ch_bam_for_calling = GATK_SPLITNCIGARREADS.out.bam
        }

        // Step 3: HaplotypeCaller
        GATK_HAPLOTYPECALLER(
            ch_bam_for_calling,
            ch_fasta.first(),
            ch_fasta_fai.first(),
            params.dbsnp ? Channel.fromPath(params.dbsnp).first() : file('NO_DBSNP'),
            params.dbsnp ? Channel.fromPath("${params.dbsnp}.tbi").first() : file('NO_DBSNP_TBI')
        )
        ch_versions = ch_versions.mix(GATK_HAPLOTYPECALLER.out.versions.first())

        // Step 4: Variant filtration
        GATK_VARIANTFILTRATION(
            GATK_HAPLOTYPECALLER.out.vcf,
            ch_fasta.first(),
            ch_fasta_fai.first()
        )
        ch_versions = ch_versions.mix(GATK_VARIANTFILTRATION.out.versions.first())

        // CHASMplus driver mutation analysis
        if (params.run_chasmplus) {
            CHASMPLUS(GATK_VARIANTFILTRATION.out.vcf, params.fasta)
        }

        // Create MAF for maftools analysis
        if (params.run_maftools) {
            CREATE_MAF(GATK_VARIANTFILTRATION.out.vcf, params.fasta, params.gtf)
        }

        // TMB estimation
        ch_metadata_file = ch_bam_bai
            .map { meta, bam, bai -> meta }
            .collect()
            .map { metas ->
                def f = file("${workDir}/metadata_for_tmb.tsv")
                f.text = "sample_id\n" + metas.collect{ it.id }.join('\n')
                return f
            }

        TMB_ESTIMATION(
            GATK_VARIANTFILTRATION.out.vcf.map{ meta, vcf, tbi -> vcf }.collect(),
            ch_metadata_file
        )
        ch_versions = ch_versions.mix(TMB_ESTIMATION.out.versions)
    }

    // ========================================
    // HLA TYPING
    // ========================================
    if (params.run_hla) {
        if (params.hla_tool == 'arcashla' || params.hla_tool == 'both') {
            ARCASHLA_EXTRACT(ch_bam_bai)
            ARCASHLA_GENOTYPE(ARCASHLA_EXTRACT.out.extracted)
            ARCASHLA_MERGE(ARCASHLA_GENOTYPE.out.genotype.map{ meta, f -> f }.collect())
            ch_versions = ch_versions.mix(ARCASHLA_EXTRACT.out.versions.first())
        }

        if (params.hla_tool == 'optitype' || params.hla_tool == 'both') {
            OPTITYPE(ch_bam_bai)
            ch_versions = ch_versions.mix(OPTITYPE.out.versions.first())
        }
    }

    emit:
    ancestry_proportions = params.run_ancestry ? ANCESTRY_INFERENCE.out.proportions : Channel.empty()
    ancestry_categories  = params.run_ancestry ? ANCESTRY_INFERENCE.out.categories : Channel.empty()
    ancestry_pca         = params.run_ancestry ? ANCESTRY_INFERENCE.out.pca : Channel.empty()
    filtered_vcf         = params.run_variant_calling ? GATK_VARIANTFILTRATION.out.vcf : Channel.empty()
    tmb_scores           = params.run_variant_calling ? TMB_ESTIMATION.out.tmb : Channel.empty()
    chasmplus_results    = params.run_variant_calling && params.run_chasmplus ? CHASMPLUS.out.results : Channel.empty()
    maf                  = params.run_variant_calling && params.run_maftools ? CREATE_MAF.out.maf : Channel.empty()
    hla_types            = params.run_hla && (params.hla_tool == 'arcashla' || params.hla_tool == 'both') ? ARCASHLA_GENOTYPE.out.genotype : Channel.empty()
    hla_merged           = params.run_hla && (params.hla_tool == 'arcashla' || params.hla_tool == 'both') ? ARCASHLA_MERGE.out.merged : Channel.empty()
    findings             = ch_findings
    versions             = ch_versions
}
