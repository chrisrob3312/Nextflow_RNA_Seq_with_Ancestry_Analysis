/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPROCESSING Subworkflow
    QC + featureCounts + normalization + batch correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FLAGSTAT            } from '../../modules/local/qc/main'
include { SAMTOOLS_IDXSTATS            } from '../../modules/local/qc/main'
include { SAMTOOLS_STATS               } from '../../modules/local/qc/main'
include { RSEQC_BAMSTAT                } from '../../modules/local/qc/main'
include { RSEQC_READDISTRIBUTION       } from '../../modules/local/qc/main'
include { RSEQC_INFEREXPERIMENT        } from '../../modules/local/qc/main'
include { PICARD_COLLECTRNASEQMETRICS  } from '../../modules/local/qc/main'
include { SUBREAD_FEATURECOUNTS        } from '../../modules/local/counting/main'
include { MERGE_COUNTS                 } from '../../modules/local/counting/main'
include { NORMALIZE_COUNTS             } from '../../modules/local/counting/main'
include { BATCH_CORRECTION             } from '../../modules/local/counting/main'
include { DETECT_GENOME_BUILD          } from '../../modules/local/genome_utils/main'
include { CROSSMAP_BAM                 } from '../../modules/local/genome_utils/main'

workflow PREPROCESSING {

    take:
    ch_bam_bai   // channel: [ val(meta), path(bam), path(bai) ]
    ch_gtf       // channel: path(gtf)
    ch_fasta     // channel: path(fasta)

    main:
    ch_versions    = Channel.empty()
    ch_qc_reports  = Channel.empty()

    // ---- Detect genome build and liftover if needed ----
    DETECT_GENOME_BUILD(ch_bam_bai)
    ch_versions = ch_versions.mix(DETECT_GENOME_BUILD.out.versions)

    // Liftover hg19 BAMs to hg38 if target genome is hg38
    if (params.target_genome_build == 'hg38') {
        ch_needs_liftover = DETECT_GENOME_BUILD.out.bam_with_build
            .branch {
                liftover: it[0].genome_build == 'hg19' || it[0].genome_build == 'GRCh37'
                pass:     true
            }

        CROSSMAP_BAM(
            ch_needs_liftover.liftover,
            ch_fasta,
            params.chain_file ? Channel.fromPath(params.chain_file).collect() : Channel.empty()
        )
        ch_versions = ch_versions.mix(CROSSMAP_BAM.out.versions)

        ch_bam_final = ch_needs_liftover.pass.mix(CROSSMAP_BAM.out.bam)
    } else {
        ch_bam_final = ch_bam_bai
    }

    // ---- QC ----
    if (params.run_qc) {
        SAMTOOLS_FLAGSTAT(ch_bam_final)
        SAMTOOLS_IDXSTATS(ch_bam_final)
        SAMTOOLS_STATS(ch_bam_final)

        ch_qc_reports = ch_qc_reports
            .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map{ meta, f -> f })
            .mix(SAMTOOLS_IDXSTATS.out.idxstats.map{ meta, f -> f })
            .mix(SAMTOOLS_STATS.out.stats.map{ meta, f -> f })

        ch_versions = ch_versions
            .mix(SAMTOOLS_FLAGSTAT.out.versions.first())
            .mix(SAMTOOLS_IDXSTATS.out.versions.first())
            .mix(SAMTOOLS_STATS.out.versions.first())

        if (!params.skip_rseqc) {
            RSEQC_BAMSTAT(ch_bam_final)
            ch_qc_reports = ch_qc_reports.mix(RSEQC_BAMSTAT.out.bam_stat.map{ meta, f -> f })
            ch_versions   = ch_versions.mix(RSEQC_BAMSTAT.out.versions.first())
        }

        if (!params.skip_picard) {
            PICARD_COLLECTRNASEQMETRICS(ch_bam_final, ch_gtf.first(), ch_fasta.first())
            ch_qc_reports = ch_qc_reports.mix(PICARD_COLLECTRNASEQMETRICS.out.metrics.map{ meta, f -> f })
            ch_versions   = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions.first())
        }
    }

    // ---- Counting ----
    if (params.run_counting) {
        SUBREAD_FEATURECOUNTS(ch_bam_final, ch_gtf.first())
        ch_qc_reports = ch_qc_reports.mix(SUBREAD_FEATURECOUNTS.out.summary.map{ meta, f -> f })
        ch_versions   = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

        // Merge individual counts into matrix
        MERGE_COUNTS(SUBREAD_FEATURECOUNTS.out.counts.map{ meta, f -> f }.collect())
        ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)

        // Write metadata to file for R scripts
        ch_metadata_file = ch_bam_final
            .map { meta, bam, bai -> meta }
            .collect()
            .map { metas ->
                def header = "sample_id\tbatch\tsex\tage\ttumor_purity\tcytomolecular_subgroup\trelapse_status\tadi_quartile\ttimepoint\tdisease_stage\tblast_percentage"
                def lines = metas.collect { m ->
                    "${m.id}\t${m.batch}\t${m.sex}\t${m.age ?: 'NA'}\t${m.tumor_purity ?: 'NA'}\t${m.cytomolecular_subgroup}\t${m.relapse_status}\t${m.adi_quartile}\t${m.timepoint}\t${m.disease_stage}\t${m.blast_percentage ?: 'NA'}"
                }
                def content = ([header] + lines).join('\n')
                def meta_file = file("${workDir}/metadata.tsv")
                meta_file.text = content
                return meta_file
            }

        // Normalize
        NORMALIZE_COUNTS(MERGE_COUNTS.out.count_matrix, ch_metadata_file)
        ch_versions = ch_versions.mix(NORMALIZE_COUNTS.out.versions)

        // Batch correction if needed
        if (params.batch_correction != 'none') {
            BATCH_CORRECTION(MERGE_COUNTS.out.count_matrix, ch_metadata_file)
            ch_versions = ch_versions.mix(BATCH_CORRECTION.out.versions)
        }
    }

    emit:
    counts_raw        = params.run_counting ? SUBREAD_FEATURECOUNTS.out.counts : Channel.empty()
    count_matrix      = params.run_counting ? MERGE_COUNTS.out.count_matrix : Channel.empty()
    counts_normalized = params.run_counting ? NORMALIZE_COUNTS.out.normalized : Channel.empty()
    qc_reports        = ch_qc_reports
    bam_final         = ch_bam_final
    versions          = ch_versions
}
