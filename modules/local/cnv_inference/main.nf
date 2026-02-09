/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNV Inference Module - InferCNV from RNA-seq expression data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process INFERCNV {
    label 'process_very_high'

    input:
    path count_matrix
    path gene_order_file
    path annotations_file  // Sample annotations (tumor vs. reference)

    output:
    path "infercnv_output/",               emit: results
    path "infercnv_output/infercnv.png",   emit: heatmap
    path "infercnv_output/infercnv.observations.txt", emit: cnv_scores
    path "versions.yml",                   emit: versions

    script:
    """
    Rscript ${projectDir}/bin/run_infercnv.R \\
        --counts ${count_matrix} \\
        --gene-order ${gene_order_file} \\
        --annotations ${annotations_file} \\
        --output-dir infercnv_output \\
        --num-threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        infercnv: \$(Rscript -e "cat(as.character(packageVersion('infercnv')))")
    END_VERSIONS
    """
}
