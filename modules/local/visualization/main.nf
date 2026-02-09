/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Visualization Module - Integrated plots and figure generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VISUALIZATION_SUMMARY {
    label 'process_medium'

    input:
    path de_results
    path pathway_results
    path immune_results
    path ancestry_results
    path wgcna_results
    path pharma_results
    path metadata

    output:
    path "figures/",             emit: figures
    path "figure_index.html",    emit: index
    path "versions.yml",         emit: versions

    script:
    """
    mkdir -p figures/{de,pathway,immune,ancestry,wgcna,pharma,overview}

    Rscript ${projectDir}/bin/create_visualizations.R \\
        --de-dir ${de_results} \\
        --pathway-dir ${pathway_results} \\
        --immune-dir ${immune_results} \\
        --ancestry-dir ${ancestry_results} \\
        --wgcna-dir ${wgcna_results} \\
        --pharma-dir ${pharma_results} \\
        --metadata ${metadata} \\
        --output-dir figures \\
        --format "png,pdf"

    # Generate HTML figure index
    python3 ${projectDir}/bin/generate_figure_index.py \\
        --figure-dir figures \\
        --output figure_index.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

process PIPELINE_LOGGER {
    label 'process_low'

    input:
    val  step_name
    val  step_description
    path input_files
    path output_summary

    output:
    path "pipeline_log_entry.md", emit: log_entry

    script:
    """
    cat <<-MDEOF > pipeline_log_entry.md
    ## \${step_name}

    **Timestamp:** \$(date -u '+%Y-%m-%d %H:%M:%S UTC')

    ### Description
    ${step_description}

    ### Input Files
    \$(ls ${input_files} 2>/dev/null | head -20 | sed 's/^/- /')

    ### Key Findings
    \$(cat ${output_summary} 2>/dev/null | head -50 || echo "See output directory for details.")

    ---

    MDEOF
    """
}

process COMPILE_PIPELINE_LOG {
    label 'process_low'

    input:
    path log_entries
    path pipeline_params

    output:
    path "pipeline_analysis_log.md", emit: log
    path "versions.yml",             emit: versions

    script:
    """
    cat <<-HEADER > pipeline_analysis_log.md
    # Cancer RNA-Seq Analysis Pipeline Log
    **Generated:** \$(date -u '+%Y-%m-%d %H:%M:%S UTC')
    **Pipeline Version:** ${workflow.manifest.version ?: 'dev'}
    **Nextflow Version:** ${nextflow.version}

    ## Pipeline Parameters
    \$(cat ${pipeline_params})

    ---

    ## Analysis Steps

    HEADER

    # Append all log entries in order
    for entry in ${log_entries}; do
        cat "\$entry" >> pipeline_analysis_log.md
    done

    cat <<-FOOTER >> pipeline_analysis_log.md

    ---

    ## Summary of Salient Findings

    *This section is auto-populated with key findings from each module.*

    ### Samples Processed
    - See QC module output for sample-level metrics

    ### Key Differential Expression Findings
    - See DE module output for significant genes per contrast

    ### Ancestry Distribution
    - See ancestry module for population assignment and proportions

    ### Notable Fusions / Structural Variants
    - See fusions module for high-confidence calls

    ### Immune Microenvironment
    - See immune deconvolution for cell type composition

    ### Drug Sensitivity Predictions
    - See pharmacogenomics module for actionable targets

    ---

    *End of pipeline log*
    FOOTER

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1)
    END_VERSIONS
    """
}
