/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sensitivity Analysis Module
    Assessments: timepoint effects, ancestry proportion thresholds,
                 covariate inclusion/exclusion, model specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SENSITIVITY_ANALYSIS {
    label 'process_high'

    input:
    path count_matrix
    path metadata
    path ancestry_proportions
    path de_results

    output:
    path "sensitivity_results/",                          emit: results
    path "sensitivity_results/timepoint_analysis.tsv",    emit: timepoint
    path "sensitivity_results/ancestry_sensitivity.tsv",  emit: ancestry
    path "sensitivity_results/covariate_impact.tsv",      emit: covariate
    path "sensitivity_results/model_comparison.tsv",      emit: model_comparison
    path "sensitivity_plots/",                            emit: plots
    path "versions.yml",                                  emit: versions

    script:
    """
    mkdir -p sensitivity_results sensitivity_plots

    Rscript ${projectDir}/bin/sensitivity_analysis.R \\
        --counts ${count_matrix} \\
        --metadata ${metadata} \\
        --ancestry ${ancestry_proportions} \\
        --de-results-dir ${de_results} \\
        --variables "${params.sensitivity_variables}" \\
        --timepoints "${params.sensitivity_timepoints}" \\
        --bootstrap-n ${params.bootstrap_iterations} \\
        --covariates "${params.de_covariates}" \\
        --output-dir sensitivity_results \\
        --plot-dir sensitivity_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}
