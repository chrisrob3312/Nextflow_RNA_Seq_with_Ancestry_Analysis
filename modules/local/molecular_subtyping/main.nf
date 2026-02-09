/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Molecular Subtyping / Cell-of-Origin Classification Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MOLECULAR_SUBTYPING {
    label 'process_medium'

    input:
    path normalized_counts
    path metadata

    output:
    path "subtyping_results/",                    emit: results
    path "subtyping_results/subtypes.tsv",        emit: subtypes
    path "subtyping_results/classifier_scores.tsv", emit: scores
    path "subtyping_plots/",                      emit: plots
    path "versions.yml",                          emit: versions

    script:
    """
    mkdir -p subtyping_results subtyping_plots

    Rscript ${projectDir}/bin/molecular_subtyping.R \\
        --expression ${normalized_counts} \\
        --metadata ${metadata} \\
        --method ${params.subtyping_method} \\
        --output-dir subtyping_results \\
        --plot-dir subtyping_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}
