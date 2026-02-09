/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pharmacogenomics Module - Drug response prediction and target ID
    Tools: oncoPredict/pRRophetic, DGIdb, CMap, DrugBank
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PHARMACOGENOMICS {
    label 'process_high'

    input:
    path normalized_counts
    path de_results
    path metadata
    path ancestry_proportions

    output:
    path "pharma_results/",                              emit: results
    path "pharma_results/drug_sensitivity_scores.tsv",   emit: drug_scores
    path "pharma_results/druggable_targets.tsv",         emit: targets
    path "pharma_results/dgidb_interactions.tsv",        emit: dgidb
    path "pharma_results/cmap_connections.tsv",          emit: cmap, optional: true
    path "pharma_results/group_comparisons/",            emit: comparisons
    path "pharma_plots/",                                emit: plots
    path "versions.yml",                                 emit: versions

    script:
    def cmap_arg = params.cmap_signatures ? "--cmap-signatures ${params.cmap_signatures}" : ''
    """
    mkdir -p pharma_results/group_comparisons pharma_plots

    Rscript ${projectDir}/bin/pharmacogenomics_analysis.R \\
        --expression ${normalized_counts} \\
        --de-results-dir ${de_results} \\
        --metadata ${metadata} \\
        --ancestry ${ancestry_proportions} \\
        --drug-db ${params.drug_response_db} \\
        --dgidb ${params.dgidb_interactions} \\
        ${cmap_arg} \\
        --output-dir pharma_results \\
        --plot-dir pharma_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        oncoPredict: \$(Rscript -e "cat(as.character(packageVersion('oncoPredict')))")
    END_VERSIONS
    """
}
