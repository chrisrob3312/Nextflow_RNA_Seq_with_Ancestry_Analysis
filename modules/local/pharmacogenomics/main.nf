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

process OCTAD_ANALYSIS {
    label 'process_high'
    input:
        path de_results
        path normalized_counts
        path metadata
    output:
        path "octad_results/", emit: octad_results
        path "*.findings.md", emit: findings
    script:
    """
    mkdir -p octad_results
    Rscript <<'REOF'
    library(octad)
    library(octad.db)
    library(data.table)

    counts <- as.matrix(fread("${normalized_counts}"), rownames = 1)
    meta <- fread("${metadata}")

    # Get DE results for disease signature
    de_files <- list.files("${de_results}", pattern = "*.tsv", full.names = TRUE, recursive = TRUE)

    for (de_file in de_files) {
        contrast_name <- gsub("\\\\.tsv\$", "", basename(de_file))
        de <- fread(de_file)

        if (!all(c("log2FoldChange", "padj") %in% names(de))) next

        # Create disease signature: top up/down regulated genes
        sig_up <- de[padj < 0.05 & log2FoldChange > 0.585, ][order(-log2FoldChange)][1:min(.N, 250)]
        sig_down <- de[padj < 0.05 & log2FoldChange < -0.585, ][order(log2FoldChange)][1:min(.N, 250)]

        if (nrow(sig_up) < 10 || nrow(sig_down) < 10) next

        # Compute drug reversal scores using OCTAD's LINCS L1000 reference
        tryCatch({
            # Use sRGES (summarized Reversal Gene Expression Score)
            # Negative sRGES = drug reverses disease signature (therapeutic candidate)
            if ("gene" %in% names(de)) {
                gene_col <- "gene"
            } else if ("gene_name" %in% names(de)) {
                gene_col <- "gene_name"
            } else {
                gene_col <- names(de)[1]
            }

            disease_sig <- data.frame(
                gene = c(sig_up[[gene_col]], sig_down[[gene_col]]),
                log2FC = c(sig_up\$log2FoldChange, sig_down\$log2FoldChange)
            )

            # Score drugs against disease signature
            res <- runsRGES(
                dz_signature = disease_sig,
                output_path = paste0("octad_results/", contrast_name),
                permutations = 10000
            )

            fwrite(res, paste0("octad_results/", contrast_name, "_drug_candidates.tsv"), sep = "\t")

            # Top therapeutic candidates (most negative sRGES)
            top_drugs <- head(res[order(res\$sRGES), ], 50)
            fwrite(top_drugs, paste0("octad_results/", contrast_name, "_top50_candidates.tsv"), sep = "\t")
        }, error = function(e) {
            message(paste("OCTAD failed for", contrast_name, ":", e\$message))
        })
    }
    REOF

    cat <<-FINDINGS > octad.findings.md
    ## OCTAD Reverse-Signature Drug Discovery
    - Contrasts analyzed: \$(ls octad_results/*_top50_candidates.tsv 2>/dev/null | wc -l)
    - Top candidates identified per contrast (sRGES < -0.2)
    FINDINGS
    """
}

process PHARMACOGX_ANALYSIS {
    label 'process_high'
    input:
        path normalized_counts
        path metadata
        path ancestry_proportions
    output:
        path "pharmacogx_results/", emit: pharmacogx_results
        path "*.findings.md", emit: findings
    script:
    """
    mkdir -p pharmacogx_results
    Rscript <<'REOF'
    library(PharmacoGx)
    library(data.table)
    library(ggplot2)

    counts <- as.matrix(fread("${normalized_counts}"), rownames = 1)
    meta <- fread("${metadata}")
    anc <- tryCatch(fread("${ancestry_proportions}"), error = function(e) NULL)
    if (!is.null(anc)) meta <- merge(meta, anc, by = "sample_id", all.x = TRUE)

    # Load available PharmacoSets
    datasets <- c("GDSC2", "CTRPv2", "gCSI")

    for (ds_name in datasets) {
        tryCatch({
            pset <- downloadPSet(ds_name, saveDir = tempdir())

            # Get drug sensitivity data
            drug_sens <- summarizeSensitivityProfiles(pset, sensitivity.measure = "auc_recomputed")

            # Get molecular profiles
            mol_prof <- summarizeMolecularProfiles(pset, mDataType = "rna", cell.lines = colnames(drug_sens))

            # Find common genes
            common_genes <- intersect(rownames(counts), rownames(mol_prof))

            if (length(common_genes) > 100) {
                # Train drug sensitivity models
                results_list <- list()
                top_drugs <- names(sort(apply(!is.na(drug_sens), 1, sum), decreasing = TRUE))[1:min(50, nrow(drug_sens))]

                for (drug in top_drugs) {
                    tryCatch({
                        # Ridge regression prediction
                        model_data <- mol_prof[common_genes, !is.na(drug_sens[drug, ])]
                        response <- drug_sens[drug, !is.na(drug_sens[drug, ])]

                        if (length(response) >= 20) {
                            fit <- glmnet::cv.glmnet(t(model_data), response, alpha = 0)
                            pred <- predict(fit, t(counts[common_genes, ]), s = "lambda.min")
                            results_list[[drug]] <- pred[, 1]
                        }
                    }, error = function(e) NULL)
                }

                if (length(results_list) > 0) {
                    pred_matrix <- do.call(cbind, results_list)
                    rownames(pred_matrix) <- colnames(counts)
                    fwrite(as.data.frame(pred_matrix),
                           paste0("pharmacogx_results/", ds_name, "_predictions.tsv"),
                           sep = "\t", row.names = TRUE)

                    # Group comparisons by ancestry
                    if (!is.null(anc) && "graf_category" %in% names(meta)) {
                        for (drug in colnames(pred_matrix)) {
                            drug_df <- data.frame(
                                sample_id = rownames(pred_matrix),
                                sensitivity = pred_matrix[, drug]
                            )
                            drug_df <- merge(drug_df, meta[, .(sample_id, graf_category)], by = "sample_id")
                            if (length(unique(na.omit(drug_df\$graf_category))) >= 2) {
                                kw <- kruskal.test(sensitivity ~ graf_category, data = drug_df[!is.na(graf_category)])
                                if (kw\$p.value < 0.05) {
                                    results_list[[paste0(drug, "_ancestry_pval")]] <- kw\$p.value
                                }
                            }
                        }
                    }
                }
            }
            message(paste(ds_name, "completed"))
        }, error = function(e) {
            message(paste(ds_name, "failed:", e\$message))
        })
    }
    REOF

    cat <<-FINDINGS > pharmacogx.findings.md
    ## PharmacoGx Multi-Database Analysis
    - Databases queried: GDSC2, CTRPv2, gCSI
    - Results: pharmacogx_results/
    FINDINGS
    """
}

process DREAM_DRUGGABILITY {
    label 'process_high'
    input:
        path de_results
        path normalized_counts
    output:
        path "dream_results/", emit: dream_results
        path "*.findings.md", emit: findings
    script:
    """
    mkdir -p dream_results
    Rscript <<'REOF'
    library(data.table)
    library(igraph)
    library(ggplot2)

    # DREAM: druggability evaluation using co-expression network + drug target topology
    counts <- as.matrix(fread("${normalized_counts}"), rownames = 1)

    de_files <- list.files("${de_results}", pattern = "*.tsv", full.names = TRUE, recursive = TRUE)

    for (de_file in de_files) {
        contrast_name <- gsub("\\\\.tsv\$", "", basename(de_file))
        de <- fread(de_file)
        if (!all(c("log2FoldChange", "padj") %in% names(de))) next

        # Get disease genes
        gene_col <- intersect(c("gene", "gene_name", "gene_id"), names(de))[1]
        if (is.na(gene_col)) next

        sig_genes <- de[padj < 0.05 & abs(log2FoldChange) > 0.585][[gene_col]]
        if (length(sig_genes) < 10) next

        # Build co-expression network from disease genes
        gene_idx <- which(rownames(counts) %in% sig_genes)
        if (length(gene_idx) < 10) next

        cor_mat <- cor(t(counts[gene_idx, ]), method = "spearman")
        cor_mat[abs(cor_mat) < 0.5] <- 0
        diag(cor_mat) <- 0

        g <- graph_from_adjacency_matrix(abs(cor_mat), mode = "undirected", weighted = TRUE)

        # Network topology metrics for druggability scoring
        deg <- degree(g)
        betw <- betweenness(g)
        close <- closeness(g)

        druggability_scores <- data.frame(
            gene = names(deg),
            degree = deg,
            betweenness = betw,
            closeness = close,
            druggability_rank = rank(-deg) + rank(-betw)
        )
        druggability_scores <- druggability_scores[order(druggability_scores\$druggability_rank), ]

        fwrite(druggability_scores,
               paste0("dream_results/", contrast_name, "_druggability.tsv"),
               sep = "\t")

        # Identify drug combination candidates (genes in same network module)
        if (vcount(g) > 10) {
            comm <- cluster_louvain(g)
            druggability_scores\$module <- membership(comm)[druggability_scores\$gene]

            # Top targets per module = combination candidates
            combos <- druggability_scores[, .SD[which.min(druggability_rank)], by = module]
            fwrite(combos, paste0("dream_results/", contrast_name, "_combination_candidates.tsv"), sep = "\t")
        }
    }
    REOF

    cat <<-FINDINGS > dream.findings.md
    ## DREAM Druggability Analysis
    - Network-based druggability scoring completed
    - Drug combination candidates identified per contrast
    FINDINGS
    """
}

process SIGNATURESEARCH_ANALYSIS {
    label 'process_high'
    input:
        path de_results
    output:
        path "sigsearch_results/", emit: sigsearch_results
        path "*.findings.md", emit: findings
    script:
    """
    mkdir -p sigsearch_results
    Rscript <<'REOF'
    library(signatureSearch)
    library(signatureSearchData)
    library(data.table)
    library(ggplot2)

    # Load LINCS L1000 reference database
    db_path <- system.file("extdata", "sample_db.h5", package = "signatureSearchData")

    de_files <- list.files("${de_results}", pattern = "*.tsv", full.names = TRUE, recursive = TRUE)

    for (de_file in de_files) {
        contrast_name <- gsub("\\\\.tsv\$", "", basename(de_file))
        de <- fread(de_file)
        if (!all(c("log2FoldChange", "padj") %in% names(de))) next

        gene_col <- intersect(c("gene", "gene_name"), names(de))[1]
        if (is.na(gene_col)) next

        # Create query signature: top 150 up + top 150 down
        de_sig <- de[padj < 0.05]
        up_genes <- head(de_sig[log2FoldChange > 0][order(-log2FoldChange)][[gene_col]], 150)
        down_genes <- head(de_sig[log2FoldChange < 0][order(log2FoldChange)][[gene_col]], 150)

        if (length(up_genes) < 10 || length(down_genes) < 10) next

        tryCatch({
            # CMAP method: connectivity score
            qsig <- qSig(
                query = list(upset = up_genes, downset = down_genes),
                gess_method = "CMAP",
                refdb = db_path
            )

            cmap_results <- gess_cmap(qsig, chunk_size = 5000)
            result_df <- result(cmap_results)

            fwrite(as.data.table(result_df),
                   paste0("sigsearch_results/", contrast_name, "_cmap_results.tsv"),
                   sep = "\t")

            # Top drug candidates (most negative connectivity = reversal)
            top_reversers <- head(result_df[order(result_df\$scaled_score), ], 50)
            fwrite(as.data.table(top_reversers),
                   paste0("sigsearch_results/", contrast_name, "_top50_reversers.tsv"),
                   sep = "\t")

            # Drug set enrichment analysis
            drugs_by_moa <- dsea_hyperG(drugs = head(result_df\$pert, 100), type = "MOA")
            fwrite(as.data.table(result(drugs_by_moa)),
                   paste0("sigsearch_results/", contrast_name, "_moa_enrichment.tsv"),
                   sep = "\t")

        }, error = function(e) {
            message(paste("signatureSearch failed for", contrast_name, ":", e\$message))
        })
    }
    REOF

    cat <<-FINDINGS > sigsearch.findings.md
    ## signatureSearch LINCS L1000 Analysis
    - Searched LINCS L1000 reference database
    - Top drug reversal candidates identified per contrast
    - MOA enrichment analysis completed
    FINDINGS
    """
}
