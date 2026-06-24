/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Differential Expression Module - DESeq2 and limma-voom
    Supports: continuous ancestry, categorical ancestry (GRAF),
              relapse status, ADI quartile, within cytomolecular subgroups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ESTIMATE_TUMOR_PURITY {
    label 'process_medium'

    input:
    path normalized_counts

    output:
    path "tumor_purity.tsv", emit: purity

    script:
    """
    Rscript <<'REOF'
    suppressPackageStartupMessages(library(estimate))

    expr <- read.delim("${normalized_counts}", row.names = 1, check.names = FALSE)

    # ESTIMATE expects a GCT-like input; write temp file
    tmp_in  <- "tmp_estimate_input.gct"
    tmp_out <- "tmp_estimate_scores.gct"

    # Write GCT format
    cat("#1.2\\n", file = tmp_in)
    cat(nrow(expr), "\\t", ncol(expr), "\\n", sep = "", file = tmp_in, append = TRUE)
    cat("Name\\tDescription\\t", paste(colnames(expr), collapse = "\\t"), "\\n",
        sep = "", file = tmp_in, append = TRUE)
    for (i in seq_len(nrow(expr))) {
        cat(rownames(expr)[i], "\\tna\\t",
            paste(expr[i, ], collapse = "\\t"), "\\n",
            sep = "", file = tmp_in, append = TRUE)
    }

    # Run ESTIMATE
    filterCommonGenes(input.f = tmp_in, output.f = "filtered.gct", id = "GeneSymbol")
    estimateScore("filtered.gct", tmp_out, platform = "${params.estimate_platform ?: 'illumina'}")

    # Parse output and extract tumor purity
    scores_raw <- read.delim(tmp_out, skip = 2, check.names = FALSE)
    rownames(scores_raw) <- scores_raw[, 1]
    scores_t <- as.data.frame(t(scores_raw[, -(1:2)]))
    colnames(scores_t) <- gsub(" ", "_", colnames(scores_t))

    # ESTIMATE tumor purity formula (Yoshihara et al.)
    if ("ESTIMATEScore" %in% colnames(scores_t)) {
        scores_t\$tumor_purity <- cos(0.6049872018 + 0.0001467884 * scores_t\$ESTIMATEScore)
        scores_t\$tumor_purity <- pmax(0, pmin(1, scores_t\$tumor_purity))
    } else {
        scores_t\$tumor_purity <- NA_real_
    }

    out <- data.frame(
        sample_id = rownames(scores_t),
        tumor_purity = round(scores_t\$tumor_purity, 4),
        stringsAsFactors = FALSE
    )
    write.table(out, "tumor_purity.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)

    cat(sprintf("Estimated tumor purity for %d samples (median: %.3f)\\n",
                sum(!is.na(out\$tumor_purity)), median(out\$tumor_purity, na.rm = TRUE)))
    REOF
    """
}

process DESEQ2_DE {
    label 'process_high'

    input:
    path count_matrix
    path metadata
    path ancestry_proportions
    path contrasts_json
    path tumor_purity

    output:
    path "deseq2_results/",    emit: results
    path "deseq2_plots/",      emit: plots
    path "deseq2_rds/",        emit: rds
    path "versions.yml",       emit: versions

    script:
    def covariates = params.de_covariates ?: 'batch,sex,age,tumor_purity'
    def subgroups  = params.cytomolecular_subgroups ?: ''
    """
    mkdir -p deseq2_results deseq2_plots deseq2_rds

    Rscript ${projectDir}/bin/run_deseq2.R \\
        --counts ${count_matrix} \\
        --metadata ${metadata} \\
        --ancestry ${ancestry_proportions} \\
        --contrasts ${contrasts_json} \\
        --tumor-purity ${tumor_purity} \\
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
    path tumor_purity

    output:
    path "limma_results/",  emit: results
    path "limma_plots/",    emit: plots
    path "limma_rds/",      emit: rds
    path "versions.yml",    emit: versions

    script:
    def covariates = params.de_covariates ?: 'batch,sex,age,tumor_purity'
    def subgroups  = params.cytomolecular_subgroups ?: ''
    """
    mkdir -p limma_results limma_plots limma_rds

    Rscript ${projectDir}/bin/run_limma_voom.R \\
        --counts ${count_matrix} \\
        --metadata ${metadata} \\
        --ancestry ${ancestry_proportions} \\
        --contrasts ${contrasts_json} \\
        --tumor-purity ${tumor_purity} \\
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
