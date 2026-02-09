/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCR/BCR Repertoire Module - TRUST4 from bulk RNA-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TRUST4 {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta

    output:
    tuple val(meta), path("${meta.id}_TRUST4_report.tsv"),          emit: report
    tuple val(meta), path("${meta.id}_TRUST4_barcode_report.tsv"),  emit: barcode_report
    tuple val(meta), path("${meta.id}_TRUST4_annot.fa"),            emit: annotations
    path "versions.yml",                                            emit: versions

    script:
    """
    run-trust4 \\
        -b ${bam} \\
        -f ${fasta} \\
        -t ${task.cpus} \\
        -o ${meta.id}_TRUST4 \\
        --ref human_IMGT+C.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trust4: \$(run-trust4 2>&1 | grep -o 'v[0-9.]*' | head -1 | sed 's/v//')
    END_VERSIONS
    """
}

process MERGE_TCR_REPORTS {
    label 'process_low'

    input:
    path trust4_reports

    output:
    path "tcr_repertoire_summary.tsv",  emit: summary
    path "tcr_diversity_metrics.tsv",   emit: diversity
    path "tcr_clonotype_tracking.tsv",  emit: clonotypes
    path "versions.yml",                emit: versions

    script:
    """
    python3 ${projectDir}/bin/merge_tcr_reports.py \\
        --input-dir . \\
        --pattern "*_TRUST4_report.tsv" \\
        --output-summary tcr_repertoire_summary.tsv \\
        --output-diversity tcr_diversity_metrics.tsv \\
        --output-clonotypes tcr_clonotype_tracking.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
