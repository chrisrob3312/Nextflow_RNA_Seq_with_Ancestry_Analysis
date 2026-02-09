/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMMUNOGENOMICS Subworkflow
    Immune deconvolution + ESTIMATE + Neoantigen prediction + TCR/BCR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IMMUNE_DECONVOLUTION } from '../../modules/local/immune_analysis/main'
include { ESTIMATE_SCORES      } from '../../modules/local/immune_analysis/main'
include { PVACSEQ              } from '../../modules/local/neoantigen/main'
include { NEOFUSE              } from '../../modules/local/neoantigen/main'
include { MERGE_NEOANTIGENS    } from '../../modules/local/neoantigen/main'
include { TRUST4               } from '../../modules/local/tcr_repertoire/main'
include { MERGE_TCR_REPORTS    } from '../../modules/local/tcr_repertoire/main'

workflow IMMUNOGENOMICS {

    take:
    ch_bam_bai           // channel: [ val(meta), path(bam), path(bai) ]
    ch_count_matrix      // channel: path(count_matrix)
    ch_normalized_counts // channel: path(normalized_counts)
    ch_variants          // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_hla_types         // channel: [ val(meta), path(hla_json) ]
    ch_fusions           // channel: [ val(meta), path(fusions) ]
    ch_metadata          // channel: val(metadata_list) collected

    main:
    ch_versions = Channel.empty()

    // Build metadata file
    ch_metadata_file = ch_metadata
        .map { metas ->
            def header = "sample_id\tbatch\tsex\tage\ttumor_purity\tcytomolecular_subgroup\trelapse_status\tadi_quartile\ttimepoint\tblast_percentage"
            def lines = metas.collect { m ->
                "${m.id}\t${m.batch}\t${m.sex}\t${m.age ?: 'NA'}\t${m.tumor_purity ?: 'NA'}\t${m.cytomolecular_subgroup}\t${m.relapse_status}\t${m.adi_quartile}\t${m.timepoint}\t${m.blast_percentage ?: 'NA'}"
            }
            def content = ([header] + lines).join('\n')
            def f = file("${workDir}/metadata_immuno.tsv")
            f.text = content
            return f
        }

    // ========================================
    // IMMUNE DECONVOLUTION
    // ========================================
    if (params.run_immune) {
        IMMUNE_DECONVOLUTION(ch_normalized_counts, ch_metadata_file)
        ESTIMATE_SCORES(ch_normalized_counts, ch_metadata_file)
        ch_versions = ch_versions
            .mix(IMMUNE_DECONVOLUTION.out.versions)
            .mix(ESTIMATE_SCORES.out.versions)
    }

    // ========================================
    // NEOANTIGEN PREDICTION
    // ========================================
    if (params.run_neoantigen && params.run_variant_calling && params.run_hla) {
        // Match VCFs with HLA types by sample
        ch_vcf_hla = ch_variants.join(ch_hla_types, by: [0])
            .map { meta, vcf, tbi, hla ->
                [ meta, vcf, tbi, hla ]
            }

        if (params.neoantigen_tool == 'pvacseq' || params.neoantigen_tool == 'both') {
            PVACSEQ(
                ch_vcf_hla.map{ meta, vcf, tbi, hla -> [ meta, vcf, tbi ] },
                ch_vcf_hla.map{ meta, vcf, tbi, hla -> [ meta, hla ] }
            )
            ch_versions = ch_versions.mix(PVACSEQ.out.versions.first())
        }

        // Fusion neoantigens via NeoFuse
        if (params.neoantigen_tool == 'neofuse' || params.neoantigen_tool == 'both') {
            if (params.run_fusions) {
                ch_fusion_hla = ch_fusions.join(ch_hla_types, by: [0])
                NEOFUSE(
                    ch_fusion_hla.map{ meta, fusions, hla -> [ meta, fusions ] },
                    ch_fusion_hla.map{ meta, fusions, hla -> [ meta, hla ] }
                )
                ch_versions = ch_versions.mix(NEOFUSE.out.versions.first())
            }
        }
    }

    // ========================================
    // TCR/BCR REPERTOIRE
    // ========================================
    if (params.run_tcr) {
        TRUST4(ch_bam_bai, Channel.fromPath(params.fasta).first())
        ch_versions = ch_versions.mix(TRUST4.out.versions.first())

        MERGE_TCR_REPORTS(TRUST4.out.report.map{ meta, f -> f }.collect())
        ch_versions = ch_versions.mix(MERGE_TCR_REPORTS.out.versions)
    }

    emit:
    immune_scores  = params.run_immune ? IMMUNE_DECONVOLUTION.out.all_scores : Channel.empty()
    estimate       = params.run_immune ? ESTIMATE_SCORES.out.scores : Channel.empty()
    neoantigens    = params.run_neoantigen && params.run_variant_calling && params.run_hla && (params.neoantigen_tool == 'pvacseq' || params.neoantigen_tool == 'both') ? PVACSEQ.out.results : Channel.empty()
    tcr_repertoire = params.run_tcr ? MERGE_TCR_REPORTS.out.summary : Channel.empty()
    versions       = ch_versions
}
