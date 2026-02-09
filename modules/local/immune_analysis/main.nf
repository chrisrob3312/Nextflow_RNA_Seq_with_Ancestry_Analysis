/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Immune Analysis Module - Deconvolution, ESTIMATE, tumor purity correction
    Tools: CIBERSORTx, xCell, MCP-counter, ESTIMATE, EPIC, TIMER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process IMMUNE_DECONVOLUTION {
    label 'process_high'

    input:
    path normalized_counts
    path metadata

    output:
    path "immune_results/",                       emit: results
    path "immune_results/deconvolution_all.tsv",  emit: all_scores
    path "immune_results/cell_fractions.tsv",     emit: cell_fractions
    path "immune_plots/",                         emit: plots
    path "versions.yml",                          emit: versions

    script:
    def cibersort_arg = params.cibersortx_token ? "--cibersortx-token ${params.cibersortx_token}" : ''
    """
    mkdir -p immune_results immune_plots

    Rscript ${projectDir}/bin/immune_deconvolution.R \\
        --expression ${normalized_counts} \\
        --metadata ${metadata} \\
        --methods "${params.immune_deconv_methods}" \\
        ${cibersort_arg} \\
        --sig-matrix ${params.cibersortx_sigmatrix} \\
        --correct-purity ${params.correct_tumor_purity} \\
        --output-dir immune_results \\
        --plot-dir immune_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        immunedeconv: \$(Rscript -e "cat(as.character(packageVersion('immunedeconv')))")
    END_VERSIONS
    """
}

process ESTIMATE_SCORES {
    label 'process_medium'

    input:
    path normalized_counts
    path metadata

    output:
    path "estimate_scores.tsv",     emit: scores
    path "estimate_purity.tsv",     emit: purity
    path "estimate_plots/",         emit: plots
    path "versions.yml",            emit: versions

    script:
    """
    mkdir -p estimate_plots

    Rscript ${projectDir}/bin/run_estimate.R \\
        --expression ${normalized_counts} \\
        --metadata ${metadata} \\
        --platform ${params.estimate_platform} \\
        --output-scores estimate_scores.tsv \\
        --output-purity estimate_purity.tsv \\
        --plot-dir estimate_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        estimate: \$(Rscript -e "cat(as.character(packageVersion('estimate')))")
    END_VERSIONS
    """
}
