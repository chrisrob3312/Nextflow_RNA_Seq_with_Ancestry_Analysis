/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Reporting Module - Integrated pipeline report generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GENERATE_REPORT {
    label 'process_medium'

    input:
    path all_results

    output:
    path "integrated_report.html",  emit: html_report
    path "summary_tables/",         emit: tables
    path "versions.yml",            emit: versions

    script:
    """
    mkdir -p summary_tables

    python3 ${projectDir}/bin/generate_report.py \\
        --results-dir . \\
        --output-html integrated_report.html \\
        --output-tables summary_tables

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
