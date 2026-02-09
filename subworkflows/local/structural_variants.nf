/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STRUCTURAL_VARIANTS Subworkflow
    Differential splicing + Gene fusions + CNV inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RMATS                     } from '../../modules/local/splicing/main'
include { LEAFCUTTER_JUNCTIONS      } from '../../modules/local/splicing/main'
include { LEAFCUTTER_CLUSTER        } from '../../modules/local/splicing/main'
include { LEAFCUTTER_DIFF_SPLICING  } from '../../modules/local/splicing/main'
include { FUSIONCATCHER             } from '../../modules/local/fusions/main'
include { ARRIBA                    } from '../../modules/local/fusions/main'
include { ARRIBA_VISUALIZATION      } from '../../modules/local/fusions/main'
include { MERGE_FUSIONS             } from '../../modules/local/fusions/main'
include { INFERCNV                  } from '../../modules/local/cnv_inference/main'

workflow STRUCTURAL_VARIANTS {

    take:
    ch_bam_bai   // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta     // channel: path(fasta)
    ch_gtf       // channel: path(gtf)
    ch_metadata  // channel: val(metadata_list) collected

    main:
    ch_versions = Channel.empty()

    // ========================================
    // DIFFERENTIAL SPLICING
    // ========================================
    if (params.run_splicing) {
        if (params.splicing_tool == 'leafcutter' || params.splicing_tool == 'both') {
            LEAFCUTTER_JUNCTIONS(ch_bam_bai)
            LEAFCUTTER_CLUSTER(
                LEAFCUTTER_JUNCTIONS.out.junctions.map{ meta, junc -> junc }.collect()
            )
            ch_versions = ch_versions.mix(LEAFCUTTER_JUNCTIONS.out.versions.first())
        }

        // rMATS requires group comparison BAM lists - handled by contrast definitions
        // The R script will generate BAM lists per contrast from metadata
    }

    // ========================================
    // FUSION DETECTION
    // ========================================
    if (params.run_fusions) {
        if (params.fusion_tool == 'arriba' || params.fusion_tool == 'both') {
            def ch_blacklist    = params.arriba_blacklist ? Channel.fromPath(params.arriba_blacklist).first() : file('NO_BLACKLIST')
            def ch_known        = params.arriba_known_fusions ? Channel.fromPath(params.arriba_known_fusions).first() : file('NO_KNOWN')

            ARRIBA(ch_bam_bai, ch_fasta.first(), ch_gtf.first(), ch_blacklist, ch_known)
            ch_versions = ch_versions.mix(ARRIBA.out.versions.first())
        }

        if (params.fusion_tool == 'fusioncatcher' || params.fusion_tool == 'both') {
            if (params.fusioncatcher_data) {
                FUSIONCATCHER(
                    ch_bam_bai,
                    Channel.fromPath(params.fusioncatcher_data).first()
                )
                ch_versions = ch_versions.mix(FUSIONCATCHER.out.versions.first())
            }
        }

        // Merge fusion results if using both tools
        if (params.fusion_tool == 'both' && params.fusioncatcher_data) {
            MERGE_FUSIONS(
                FUSIONCATCHER.out.fusions.map{ meta, f -> f }.collect(),
                ARRIBA.out.fusions.map{ meta, f -> f }.collect()
            )
            ch_versions = ch_versions.mix(MERGE_FUSIONS.out.versions)
        }
    }

    // ========================================
    // CNV INFERENCE FROM RNA-SEQ
    // ========================================
    if (params.run_cnv_inference) {
        if (params.cnv_gene_order_file) {
            // InferCNV needs a gene order file and reference/observation annotations
            INFERCNV(
                Channel.empty(),  // Count matrix fed separately
                Channel.fromPath(params.cnv_gene_order_file).first(),
                Channel.empty()   // Annotations generated from metadata
            )
            ch_versions = ch_versions.mix(INFERCNV.out.versions)
        }
    }

    emit:
    splicing_results = params.run_splicing && (params.splicing_tool == 'leafcutter' || params.splicing_tool == 'both') ? LEAFCUTTER_CLUSTER.out.counts : Channel.empty()
    fusion_results   = params.run_fusions && (params.fusion_tool == 'arriba' || params.fusion_tool == 'both') ? ARRIBA.out.fusions : Channel.empty()
    cnv_results      = params.run_cnv_inference && params.cnv_gene_order_file ? INFERCNV.out.results : Channel.empty()
    versions         = ch_versions
}
