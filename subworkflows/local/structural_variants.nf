/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STRUCTURAL_VARIANTS Subworkflow
    Differential splicing (rMATS + Leafcutter + SplAdder + Bisbee)
    Gene fusions + CNV inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RMATS                     } from '../../modules/local/splicing/main'
include { LEAFCUTTER_JUNCTIONS      } from '../../modules/local/splicing/main'
include { LEAFCUTTER_CLUSTER        } from '../../modules/local/splicing/main'
include { LEAFCUTTER_DIFF_SPLICING  } from '../../modules/local/splicing/main'
include { SPLADDER_BUILD            } from '../../modules/local/splicing/main'
include { SPLADDER_TEST             } from '../../modules/local/splicing/main'
include { BISBEE_PREP               } from '../../modules/local/splicing/main'
include { BISBEE_DIFF               } from '../../modules/local/splicing/main'
include { BISBEE_PROT               } from '../../modules/local/splicing/main'
include { BISBEE_OUTLIER            } from '../../modules/local/splicing/main'
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
    ch_versions  = Channel.empty()
    ch_findings  = Channel.empty()

    // ========================================
    // DIFFERENTIAL SPLICING
    // ========================================
    if (params.run_splicing) {
        // --- Leafcutter ---
        if (params.splicing_tool in ['leafcutter', 'both', 'all']) {
            LEAFCUTTER_JUNCTIONS(ch_bam_bai)
            LEAFCUTTER_CLUSTER(
                LEAFCUTTER_JUNCTIONS.out.junctions.map{ meta, junc -> junc }.collect()
            )
            ch_versions = ch_versions.mix(LEAFCUTTER_JUNCTIONS.out.versions.first())
        }

        // --- rMATS ---
        // rMATS requires group comparison BAM lists - handled by contrast definitions
        // The contrast-specific BAM lists are generated from metadata

        // --- SplAdder + Bisbee ---
        if (params.splicing_tool in ['spladder', 'bisbee', 'all']) {
            // Collect all BAMs and BAIs for SplAdder graph building
            ch_all_bams = ch_bam_bai.map{ meta, bam, bai -> bam }.collect()
            ch_all_bais = ch_bam_bai.map{ meta, bam, bai -> bai }.collect()

            SPLADDER_BUILD(
                ch_all_bams,
                ch_all_bais,
                ch_gtf.first()
            )
            ch_versions = ch_versions.mix(SPLADDER_BUILD.out.versions)

            // Prepare Bisbee input from SplAdder counts
            BISBEE_PREP(SPLADDER_BUILD.out.graph)

            // Bisbee outlier detection (runs on all samples, no contrast needed)
            BISBEE_OUTLIER(BISBEE_PREP.out.prepped)
            ch_findings = ch_findings.mix(BISBEE_OUTLIER.out.findings)
            ch_versions = ch_versions.mix(BISBEE_OUTLIER.out.versions)
        }
    }

    // ========================================
    // FUSION DETECTION
    // ========================================
    if (params.run_fusions) {
        if (params.fusion_tool == 'arriba' || params.fusion_tool == 'both') {
            def ch_blacklist = params.arriba_blacklist ? Channel.fromPath(params.arriba_blacklist).first() : file('NO_BLACKLIST')
            def ch_known     = params.arriba_known_fusions ? Channel.fromPath(params.arriba_known_fusions).first() : file('NO_KNOWN')

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
            INFERCNV(
                Channel.empty(),  // Count matrix fed separately
                Channel.fromPath(params.cnv_gene_order_file).first(),
                Channel.empty()   // Annotations generated from metadata
            )
            ch_versions = ch_versions.mix(INFERCNV.out.versions)
        }
    }

    emit:
    splicing_results = params.run_splicing && (params.splicing_tool in ['leafcutter', 'both', 'all']) ? LEAFCUTTER_CLUSTER.out.counts : Channel.empty()
    spladder_results = params.run_splicing && (params.splicing_tool in ['spladder', 'bisbee', 'all']) ? SPLADDER_BUILD.out.graph : Channel.empty()
    bisbee_input     = params.run_splicing && (params.splicing_tool in ['spladder', 'bisbee', 'all']) ? BISBEE_PREP.out.prepped : Channel.empty()
    fusion_results   = params.run_fusions && (params.fusion_tool == 'arriba' || params.fusion_tool == 'both') ? ARRIBA.out.fusions : Channel.empty()
    cnv_results      = params.run_cnv_inference && params.cnv_gene_order_file ? INFERCNV.out.results : Channel.empty()
    findings         = ch_findings
    versions         = ch_versions
}
