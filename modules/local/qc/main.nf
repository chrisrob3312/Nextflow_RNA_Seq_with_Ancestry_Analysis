/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC Module - BAM quality control and metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SAMTOOLS_FLAGSTAT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    path "versions.yml",                 emit: versions

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.id}.flagstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}

process SAMTOOLS_IDXSTATS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.idxstats"), emit: idxstats
    path "versions.yml",                 emit: versions

    script:
    """
    samtools idxstats ${bam} > ${meta.id}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}

process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path "versions.yml",              emit: versions

    script:
    """
    samtools stats -@ ${task.cpus} ${bam} > ${meta.id}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}

process RSEQC_BAMSTAT {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam_stat.txt"), emit: bam_stat
    path "versions.yml",                     emit: versions

    script:
    """
    bam_stat.py -i ${bam} > ${meta.id}.bam_stat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(bam_stat.py --version 2>&1 | sed 's/.*bam_stat.py //')
    END_VERSIONS
    """
}

process RSEQC_READDISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path  gene_bed

    output:
    tuple val(meta), path("*.read_distribution.txt"), emit: read_dist
    path "versions.yml",                              emit: versions

    script:
    """
    read_distribution.py -i ${bam} -r ${gene_bed} > ${meta.id}.read_distribution.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_distribution.py --version 2>&1 | sed 's/.*read_distribution.py //')
    END_VERSIONS
    """
}

process RSEQC_INFEREXPERIMENT {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path  gene_bed

    output:
    tuple val(meta), path("*.infer_experiment.txt"), emit: infer_exp
    path "versions.yml",                             emit: versions

    script:
    """
    infer_experiment.py -i ${bam} -r ${gene_bed} > ${meta.id}.infer_experiment.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version 2>&1 | sed 's/.*infer_experiment.py //')
    END_VERSIONS
    """
}

process PICARD_COLLECTRNASEQMETRICS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  gtf

    output:
    tuple val(meta), path("*.rna_metrics"), emit: metrics
    path "versions.yml",                    emit: versions

    script:
    """
    # Generate refFlat from GTF
    gtfToGenePred -genePredExt ${gtf} /dev/stdout | \\
        awk 'BEGIN{OFS="\\t"}{print \$12,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' > refFlat.txt || true

    picard CollectRnaSeqMetrics \\
        -I ${bam} \\
        -O ${meta.id}.rna_metrics \\
        -REF_FLAT refFlat.txt \\
        -STRAND SECOND_READ_TRANSCRIPTION_STRAND \\
        -R ${fasta} \\
        --VALIDATION_STRINGENCY LENIENT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectRnaSeqMetrics --version 2>&1 | sed 's/.*Version://' || echo 'unknown')
    END_VERSIONS
    """
}

process MULTIQC {
    label 'process_medium'

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data",               emit: data
    path "versions.yml",         emit: versions

    script:
    def args = params.multiqc_config ? "--config ${params.multiqc_config}" : ''
    """
    multiqc . \\
        ${args} \\
        --force \\
        --title "Cancer RNA-Seq Pipeline QC Report"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/.*version //')
    END_VERSIONS
    """
}
