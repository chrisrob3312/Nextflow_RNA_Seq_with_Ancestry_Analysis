/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINICAL_ANALYSIS Subworkflow
    Sensitivity analysis + Pharmacogenomics + Molecular subtyping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SENSITIVITY_ANALYSIS } from '../../modules/local/sensitivity/main'
include { PHARMACOGENOMICS         } from '../../modules/local/pharmacogenomics/main'
include { OCTAD_ANALYSIS           } from '../../modules/local/pharmacogenomics/main'
include { PHARMACOGX_ANALYSIS      } from '../../modules/local/pharmacogenomics/main'
include { DREAM_DRUGGABILITY       } from '../../modules/local/pharmacogenomics/main'
include { SIGNATURESEARCH_ANALYSIS } from '../../modules/local/pharmacogenomics/main'
include { MOLECULAR_SUBTYPING     } from '../../modules/local/molecular_subtyping/main'

workflow CLINICAL_ANALYSIS {

    take:
    ch_count_matrix       // channel: path(count_matrix)
    ch_normalized_counts  // channel: path(normalized_counts)
    ch_de_results         // channel: path(de_results_dir)
    ch_ancestry           // channel: path(ancestry_proportions)
    ch_immune             // channel: path(immune_scores)
    ch_metadata           // channel: val(metadata_list) collected

    main:
    ch_versions = Channel.empty()

    // Build metadata file
    ch_metadata_file = ch_metadata
        .map { metas ->
            def header = "sample_id\tbatch\tsex\tage\ttumor_purity\tcytomolecular_subgroup\trelapse_status\tadi_quartile\ttimepoint\tdisease_stage\tblast_percentage"
            def lines = metas.collect { m ->
                "${m.id}\t${m.batch}\t${m.sex}\t${m.age ?: 'NA'}\t${m.tumor_purity ?: 'NA'}\t${m.cytomolecular_subgroup}\t${m.relapse_status}\t${m.adi_quartile}\t${m.timepoint}\t${m.disease_stage}\t${m.blast_percentage ?: 'NA'}"
            }
            def content = ([header] + lines).join('\n')
            def f = file("${workDir}/metadata_clinical.tsv")
            f.text = content
            return f
        }

    // ========================================
    // SENSITIVITY ANALYSIS
    // ========================================
    if (params.run_sensitivity) {
        SENSITIVITY_ANALYSIS(
            ch_count_matrix,
            ch_metadata_file,
            ch_ancestry,
            ch_de_results
        )
        ch_versions = ch_versions.mix(SENSITIVITY_ANALYSIS.out.versions)
    }

    // ========================================
    // PHARMACOGENOMICS
    // ========================================
    if (params.run_pharmacogenomics) {
        PHARMACOGENOMICS(
            ch_normalized_counts,
            ch_de_results,
            ch_metadata_file,
            ch_ancestry
        )
        ch_versions = ch_versions.mix(PHARMACOGENOMICS.out.versions)
    }

    // ========================================
    // EXTENDED PHARMACOGENOMICS
    // ========================================
    if (params.run_octad) {
        OCTAD_ANALYSIS(ch_de_results, ch_normalized_counts, ch_metadata_file)
        ch_versions = ch_versions.mix(OCTAD_ANALYSIS.out.versions)
    }

    if (params.run_pharmacogx) {
        PHARMACOGX_ANALYSIS(ch_normalized_counts, ch_metadata_file, ch_ancestry)
        ch_versions = ch_versions.mix(PHARMACOGX_ANALYSIS.out.versions)
    }

    if (params.run_dream) {
        DREAM_DRUGGABILITY(ch_de_results, ch_normalized_counts)
        ch_versions = ch_versions.mix(DREAM_DRUGGABILITY.out.versions)
    }

    if (params.run_signaturesearch) {
        SIGNATURESEARCH_ANALYSIS(ch_de_results)
        ch_versions = ch_versions.mix(SIGNATURESEARCH_ANALYSIS.out.versions)
    }

    // ========================================
    // MOLECULAR SUBTYPING
    // ========================================
    if (params.run_molecular_subtyping) {
        MOLECULAR_SUBTYPING(
            ch_normalized_counts,
            ch_metadata_file
        )
        ch_versions = ch_versions.mix(MOLECULAR_SUBTYPING.out.versions)
    }

    emit:
    sensitivity_results = params.run_sensitivity ? SENSITIVITY_ANALYSIS.out.results : Channel.empty()
    pharma_results      = params.run_pharmacogenomics ? PHARMACOGENOMICS.out.results : Channel.empty()
    subtype_results     = params.run_molecular_subtyping ? MOLECULAR_SUBTYPING.out.results : Channel.empty()
    versions            = ch_versions
}
