/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fusion Detection Module - FusionCatcher and Arriba
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FUSIONCATCHER {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fusioncatcher_data

    output:
    tuple val(meta), path("${meta.id}_fusioncatcher/"),                          emit: results
    tuple val(meta), path("${meta.id}_fusioncatcher/final-list_candidate-fusion-genes.txt"), emit: fusions
    path "versions.yml",                                                          emit: versions

    script:
    """
    mkdir -p ${meta.id}_fusioncatcher

    # Convert BAM to FASTQ for FusionCatcher
    samtools sort -n -@ ${task.cpus} ${bam} | \\
        samtools fastq -@ ${task.cpus} \\
            -1 ${meta.id}_R1.fastq.gz \\
            -2 ${meta.id}_R2.fastq.gz \\
            -0 /dev/null \\
            -s /dev/null \\
            -

    fusioncatcher.py \\
        -d ${fusioncatcher_data} \\
        -i ${meta.id}_R1.fastq.gz,${meta.id}_R2.fastq.gz \\
        -o ${meta.id}_fusioncatcher \\
        -p ${task.cpus} \\
        --skip-blat

    # Cleanup intermediate FASTQs
    rm -f ${meta.id}_R1.fastq.gz ${meta.id}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(fusioncatcher.py --version 2>&1 | head -1 | sed 's/.*version //')
    END_VERSIONS
    """
}

process ARRIBA {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  gtf
    path  blacklist
    path  known_fusions

    output:
    tuple val(meta), path("${meta.id}.arriba.fusions.tsv"),        emit: fusions
    tuple val(meta), path("${meta.id}.arriba.fusions.discarded.tsv"), emit: discarded
    path "versions.yml",                                            emit: versions

    script:
    def blacklist_arg = blacklist ? "-b ${blacklist}" : ''
    def known_arg     = known_fusions ? "-k ${known_fusions}" : ''
    """
    arriba \\
        -x ${bam} \\
        -a ${fasta} \\
        -g ${gtf} \\
        ${blacklist_arg} \\
        ${known_arg} \\
        -o ${meta.id}.arriba.fusions.tsv \\
        -O ${meta.id}.arriba.fusions.discarded.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h 2>&1 | grep 'Version' | sed 's/.*: //')
    END_VERSIONS
    """
}

process ARRIBA_VISUALIZATION {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(fusions)
    path  gtf
    path  protein_domains

    output:
    tuple val(meta), path("${meta.id}.arriba.pdf"), emit: pdf
    path "versions.yml",                            emit: versions

    script:
    """
    draw_fusions.R \\
        --fusions=${fusions} \\
        --alignments=${bam} \\
        --annotation=${gtf} \\
        --proteinDomains=${protein_domains} \\
        --output=${meta.id}.arriba.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h 2>&1 | grep 'Version' | sed 's/.*: //')
    END_VERSIONS
    """
}

process MERGE_FUSIONS {
    label 'process_low'

    input:
    path fusioncatcher_results
    path arriba_results

    output:
    path "merged_fusions.tsv",         emit: merged
    path "high_confidence_fusions.tsv", emit: high_confidence
    path "fusion_summary.tsv",         emit: summary
    path "versions.yml",               emit: versions

    script:
    """
    python3 ${projectDir}/bin/parse_fusions.py \\
        --fusioncatcher-dir . \\
        --arriba-dir . \\
        --min-reads ${params.min_fusion_reads} \\
        --output-merged merged_fusions.tsv \\
        --output-hc high_confidence_fusions.tsv \\
        --output-summary fusion_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
