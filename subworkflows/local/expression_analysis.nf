/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXPRESSION_ANALYSIS Subworkflow
    DE (DESeq2/limma) + WGCNA + Pathway Enrichment + GSVA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ESTIMATE_TUMOR_PURITY } from '../../modules/local/differential_expression/main'
include { DESEQ2_DE             } from '../../modules/local/differential_expression/main'
include { LIMMA_VOOM_DE         } from '../../modules/local/differential_expression/main'
include { WGCNA_ANALYSIS     } from '../../modules/local/wgcna/main'
include { PATHWAY_ENRICHMENT } from '../../modules/local/pathway_enrichment/main'
include { GSVA_ANALYSIS      } from '../../modules/local/pathway_enrichment/main'

workflow EXPRESSION_ANALYSIS {

    take:
    ch_count_matrix       // channel: path(raw_count_matrix.tsv)
    ch_normalized_counts  // channel: path(normalized_counts.tsv)
    ch_metadata           // channel: val(metadata_list) collected
    ch_ancestry           // channel: path(ancestry_proportions.tsv)

    main:
    ch_versions = Channel.empty()

    // Build metadata file from collected metadata
    ch_metadata_file = ch_metadata
        .map { metas ->
            def header = "sample_id\tbatch\tsex\tage\ttumor_purity\tcytomolecular_subgroup\trelapse_status\tadi_quartile\ttimepoint\tdisease_stage\tblast_percentage"
            def lines = metas.collect { m ->
                "${m.id}\t${m.batch}\t${m.sex}\t${m.age ?: 'NA'}\t${m.tumor_purity ?: 'NA'}\t${m.cytomolecular_subgroup}\t${m.relapse_status}\t${m.adi_quartile}\t${m.timepoint}\t${m.disease_stage}\t${m.blast_percentage ?: 'NA'}"
            }
            def content = ([header] + lines).join('\n')
            def f = file("${workDir}/metadata_expr.tsv")
            f.text = content
            return f
        }

    // Contrasts JSON
    ch_contrasts = params.de_contrasts
        ? Channel.fromPath(params.de_contrasts).first()
        : Channel.fromPath("${projectDir}/assets/default_contrasts.json").first()

    // ========================================
    // ESTIMATE TUMOR PURITY (data-derived, runs before DE)
    // ========================================
    ESTIMATE_TUMOR_PURITY(ch_normalized_counts)

    // ========================================
    // DIFFERENTIAL EXPRESSION
    // ========================================
    if (params.run_de) {
        if (params.de_tool == 'deseq2' || params.de_tool == 'both') {
            DESEQ2_DE(
                ch_count_matrix,
                ch_metadata_file,
                ch_ancestry,
                ch_contrasts,
                ESTIMATE_TUMOR_PURITY.out.purity
            )
            ch_versions = ch_versions.mix(DESEQ2_DE.out.versions)
        }

        if (params.de_tool == 'limma' || params.de_tool == 'both') {
            LIMMA_VOOM_DE(
                ch_count_matrix,
                ch_metadata_file,
                ch_ancestry,
                ch_contrasts,
                ESTIMATE_TUMOR_PURITY.out.purity
            )
            ch_versions = ch_versions.mix(LIMMA_VOOM_DE.out.versions)
        }
    }

    // ========================================
    // WGCNA
    // ========================================
    if (params.run_wgcna) {
        WGCNA_ANALYSIS(
            ch_normalized_counts,
            ch_metadata_file,
            ch_ancestry
        )
        ch_versions = ch_versions.mix(WGCNA_ANALYSIS.out.versions)
    }

    // ========================================
    // PATHWAY ENRICHMENT
    // ========================================
    if (params.run_pathway) {
        // Collect DE results directory
        ch_de_for_pathway = Channel.empty()
        if (params.de_tool == 'deseq2' || params.de_tool == 'both') {
            ch_de_for_pathway = DESEQ2_DE.out.results
        } else if (params.de_tool == 'limma') {
            ch_de_for_pathway = LIMMA_VOOM_DE.out.results
        }

        PATHWAY_ENRICHMENT(ch_de_for_pathway, ch_count_matrix)
        ch_versions = ch_versions.mix(PATHWAY_ENRICHMENT.out.versions)

        // GSVA
        if (params.run_gsva) {
            GSVA_ANALYSIS(ch_normalized_counts, ch_metadata_file)
            ch_versions = ch_versions.mix(GSVA_ANALYSIS.out.versions)
        }
    }

    emit:
    tumor_purity    = ESTIMATE_TUMOR_PURITY.out.purity
    de_results      = params.run_de && (params.de_tool == 'deseq2' || params.de_tool == 'both') ? DESEQ2_DE.out.results : (params.run_de && params.de_tool == 'limma' ? LIMMA_VOOM_DE.out.results : Channel.empty())
    wgcna_modules   = params.run_wgcna ? WGCNA_ANALYSIS.out.results : Channel.empty()
    pathway_results = params.run_pathway ? PATHWAY_ENRICHMENT.out.results : Channel.empty()
    gsva_scores     = params.run_pathway && params.run_gsva ? GSVA_ANALYSIS.out.scores : Channel.empty()
    versions        = ch_versions
}
