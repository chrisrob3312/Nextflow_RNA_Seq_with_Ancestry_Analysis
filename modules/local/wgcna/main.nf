/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WGCNA Module - Weighted Gene Co-expression Network Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process WGCNA_ANALYSIS {
    label 'process_very_high'

    input:
    path normalized_counts
    path metadata
    path ancestry_proportions

    output:
    path "wgcna_results/",                          emit: results
    path "wgcna_results/module_eigengenes.tsv",     emit: eigengenes
    path "wgcna_results/module_membership.tsv",     emit: membership
    path "wgcna_results/module_trait_cor.tsv",       emit: trait_correlations
    path "wgcna_results/hub_genes.tsv",             emit: hub_genes
    path "wgcna_plots/",                            emit: plots
    path "versions.yml",                            emit: versions

    script:
    """
    mkdir -p wgcna_results wgcna_plots

    Rscript ${projectDir}/bin/run_wgcna.R \\
        --expression ${normalized_counts} \\
        --metadata ${metadata} \\
        --ancestry ${ancestry_proportions} \\
        --min-module-size ${params.wgcna_min_module_size} \\
        --merge-cut-height ${params.wgcna_merge_cut_height} \\
        --network-type ${params.wgcna_network_type} \\
        --tom-type ${params.wgcna_tom_type} \\
        ${params.wgcna_soft_power ? "--soft-power ${params.wgcna_soft_power}" : ''} \\
        --output-dir wgcna_results \\
        --plot-dir wgcna_plots \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        WGCNA: \$(Rscript -e "cat(as.character(packageVersion('WGCNA')))")
    END_VERSIONS
    """
}
