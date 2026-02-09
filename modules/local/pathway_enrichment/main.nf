/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pathway Enrichment Module - ORA, GSEA, GSVA
    Databases: GO, KEGG, Reactome, Hallmark, ImmuneSigDB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PATHWAY_ENRICHMENT {
    label 'process_high'

    input:
    path de_results  // Directory with DE result tables
    path count_matrix

    output:
    path "pathway_results/",              emit: results
    path "pathway_results/ora/",          emit: ora
    path "pathway_results/gsea/",         emit: gsea
    path "pathway_plots/",               emit: plots
    path "versions.yml",                  emit: versions

    script:
    """
    mkdir -p pathway_results/ora pathway_results/gsea pathway_plots

    Rscript ${projectDir}/bin/pathway_enrichment.R \\
        --de-results-dir ${de_results} \\
        --counts ${count_matrix} \\
        --databases "${params.pathway_databases}" \\
        --species "${params.msigdb_species}" \\
        --min-size ${params.gsea_min_size} \\
        --max-size ${params.gsea_max_size} \\
        --padj-threshold ${params.padj_threshold} \\
        --output-dir pathway_results \\
        --plot-dir pathway_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        clusterProfiler: \$(Rscript -e "cat(as.character(packageVersion('clusterProfiler')))")
        fgsea: \$(Rscript -e "cat(as.character(packageVersion('fgsea')))")
        msigdbr: \$(Rscript -e "cat(as.character(packageVersion('msigdbr')))")
    END_VERSIONS
    """
}

process GSVA_ANALYSIS {
    label 'process_high'

    input:
    path normalized_counts
    path metadata

    output:
    path "gsva_scores.tsv",       emit: scores
    path "gsva_results/",         emit: results
    path "gsva_plots/",           emit: plots
    path "versions.yml",          emit: versions

    script:
    """
    mkdir -p gsva_results gsva_plots

    Rscript ${projectDir}/bin/gsva_analysis.R \\
        --expression ${normalized_counts} \\
        --metadata ${metadata} \\
        --databases "${params.pathway_databases}" \\
        --species "${params.msigdb_species}" \\
        --min-size ${params.gsea_min_size} \\
        --max-size ${params.gsea_max_size} \\
        --output-scores gsva_scores.tsv \\
        --output-dir gsva_results \\
        --plot-dir gsva_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        GSVA: \$(Rscript -e "cat(as.character(packageVersion('GSVA')))")
    END_VERSIONS
    """
}
