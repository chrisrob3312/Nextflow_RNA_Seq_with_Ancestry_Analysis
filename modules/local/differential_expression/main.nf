/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Differential Expression Module - DESeq2 and limma-voom
    Supports: continuous ancestry, categorical ancestry (GRAF),
              relapse status, ADI quartile, within cytomolecular subgroups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DESEQ2_DE {
    label 'process_high'

    input:
    path count_matrix
    path metadata
    path ancestry_proportions
    path contrasts_json

    output:
    path "deseq2_results/",    emit: results
    path "deseq2_plots/",      emit: plots
    path "deseq2_rds/",        emit: rds
    path "versions.yml",       emit: versions

    script:
    def covariates = params.de_covariates ?: 'batch,sex,age,blast_percentage'
    def subgroups  = params.cytomolecular_subgroups ?: ''
    """
    mkdir -p deseq2_results deseq2_plots deseq2_rds

    Rscript ${projectDir}/bin/run_deseq2.R \\
        --counts ${count_matrix} \\
        --metadata ${metadata} \\
        --ancestry ${ancestry_proportions} \\
        --contrasts ${contrasts_json} \\
        --covariates "${covariates}" \\
        --cytomolecular-subgroups "${subgroups}" \\
        --padj-threshold ${params.padj_threshold} \\
        --lfc-threshold ${params.lfc_threshold} \\
        --output-dir deseq2_results \\
        --plot-dir deseq2_plots \\
        --rds-dir deseq2_rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        DESeq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}

process LIMMA_VOOM_DE {
    label 'process_high'

    input:
    path count_matrix
    path metadata
    path ancestry_proportions
    path contrasts_json

    output:
    path "limma_results/",  emit: results
    path "limma_plots/",    emit: plots
    path "limma_rds/",      emit: rds
    path "versions.yml",    emit: versions

    script:
    def covariates = params.de_covariates ?: 'batch,sex,age,blast_percentage'
    def subgroups  = params.cytomolecular_subgroups ?: ''
    """
    mkdir -p limma_results limma_plots limma_rds

    Rscript ${projectDir}/bin/run_limma_voom.R \\
        --counts ${count_matrix} \\
        --metadata ${metadata} \\
        --ancestry ${ancestry_proportions} \\
        --contrasts ${contrasts_json} \\
        --covariates "${covariates}" \\
        --cytomolecular-subgroups "${subgroups}" \\
        --padj-threshold ${params.padj_threshold} \\
        --lfc-threshold ${params.lfc_threshold} \\
        --output-dir limma_results \\
        --plot-dir limma_plots \\
        --rds-dir limma_rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        limma: \$(Rscript -e "cat(as.character(packageVersion('limma')))")
        edgeR: \$(Rscript -e "cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """
}
