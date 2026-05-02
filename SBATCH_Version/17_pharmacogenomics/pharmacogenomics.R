#!/usr/bin/env Rscript
# ============================================================================
# 17 - Pharmacogenomics Analysis
# ============================================================================
# Predicts drug sensitivity using oncoPredict (GDSC/CCLE training data),
# queries DGIdb API for druggable targets from DE gene lists, and compares
# drug sensitivity scores between ancestry/clinical groups (Wilcoxon tests).
#
# Required packages: oncoPredict, httr, jsonlite, ggplot2, tidyverse, optparse
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(oncoPredict)
    library(tidyverse)
    library(ggplot2)
    library(httr)
    library(jsonlite)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
    make_option("--normalized_counts", type = "character",
                help = "Path to normalized count matrix (genes x samples TSV)"),
    make_option("--de_results_dir", type = "character",
                help = "Directory containing DE results"),
    make_option("--metadata", type = "character",
                help = "Path to sample metadata CSV"),
    make_option("--ancestry_proportions", type = "character", default = "none",
                help = "Path to ancestry proportions TSV [default: none]"),
    make_option("--output_dir", type = "character",
                help = "Output directory"),
    make_option("--drug_db", type = "character", default = "GDSC",
                help = "Drug response database: GDSC or CCLE [default: %default]"),
    make_option("--threads", type = "integer", default = 8,
                help = "Number of threads [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$normalized_counts) || is.null(opt$de_results_dir) ||
    is.null(opt$metadata) || is.null(opt$output_dir)) {
    stop("Required: --normalized_counts, --de_results_dir, --metadata, --output_dir")
}

cat("=== Pharmacogenomics Analysis ===\n")
cat("Normalized counts:", opt$normalized_counts, "\n")
cat("DE results dir:", opt$de_results_dir, "\n")
cat("Metadata:", opt$metadata, "\n")
cat("Ancestry proportions:", opt$ancestry_proportions, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Drug database:", opt$drug_db, "\n")
cat("Threads:", opt$threads, "\n\n")

# ============================================================================
# Load data
# ============================================================================
cat("Loading data...\n")

counts <- read.table(opt$normalized_counts, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE)
metadata <- read.csv(opt$metadata, header = TRUE, stringsAsFactors = FALSE)

cat(sprintf("  Count matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  Metadata: %d samples\n", nrow(metadata)))

# Load ancestry if available
ancestry <- NULL
if (opt$ancestry_proportions != "none" && file.exists(opt$ancestry_proportions)) {
    ancestry <- read.table(opt$ancestry_proportions, header = TRUE, sep = "\t",
                           row.names = 1, check.names = FALSE)
    cat(sprintf("  Ancestry: %d samples x %d components\n", nrow(ancestry), ncol(ancestry)))
}

# Load DE results
de_files <- list.files(opt$de_results_dir, pattern = "\\.tsv$",
                       full.names = TRUE, recursive = TRUE)
cat(sprintf("  Found %d DE result files\n", length(de_files)))

# Collect significant DE genes
all_de_genes <- character(0)
for (de_file in de_files) {
    tryCatch({
        de_df <- read.table(de_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        if ("padj" %in% colnames(de_df) && "log2FoldChange" %in% colnames(de_df)) {
            sig <- de_df %>%
                filter(padj < 0.05, abs(log2FoldChange) > 0.585) %>%
                pull(1)
            all_de_genes <- unique(c(all_de_genes, sig))
        }
    }, error = function(e) NULL)
}
cat(sprintf("  Total significant DE genes: %d\n", length(all_de_genes)))

# Align metadata with counts
sample_ids <- colnames(counts)
metadata <- metadata %>% filter(metadata[[1]] %in% sample_ids)
rownames(metadata) <- metadata[[1]]

# ============================================================================
# 1. Drug Sensitivity Prediction (oncoPredict)
# ============================================================================
cat("\n--- 1. Drug Sensitivity Prediction ---\n")

expr_matrix <- as.matrix(counts)

tryCatch({
    # oncoPredict uses GDSC2 or CCLE training data
    # The training data should be available via oncoPredict package
    if (toupper(opt$drug_db) == "GDSC") {
        cat("  Using GDSC2 training data...\n")
        # Load GDSC training expression and drug response
        # oncoPredict provides these as built-in data
        data("GDSC2_Expr", package = "oncoPredict", envir = environment())
        data("GDSC2_Res", package = "oncoPredict", envir = environment())
        train_expr <- GDSC2_Expr
        train_resp <- GDSC2_Res
    } else {
        cat("  Using CCLE training data...\n")
        data("CCLE_Expr", package = "oncoPredict", envir = environment())
        data("CCLE_Res", package = "oncoPredict", envir = environment())
        train_expr <- CCLE_Expr
        train_resp <- CCLE_Res
    }

    # Run calcPhenotype for drug sensitivity prediction
    cat("  Running calcPhenotype (ridge regression)...\n")
    sensitivity_dir <- file.path(opt$output_dir, "drug_sensitivity")
    dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

    calcPhenotype(
        trainingExprData = train_expr,
        trainingPtype = train_resp,
        testExprData = expr_matrix,
        batchCorrect = "eb",
        powerTransformPhenotype = TRUE,
        removeLowVaryingGenes = 0.2,
        minNumSamples = 10,
        printOutput = TRUE,
        removeLowVaringGenesFrom = "rawData",
        report_pc = FALSE,
        cc = TRUE,
        rsq = FALSE,
        pcr = FALSE,
        selection = -1
    )

    # Read predicted drug sensitivity (calcPhenotype writes to working dir)
    drug_predictions_file <- "calcPhenotype_Output/DrugPredictions.csv"
    if (file.exists(drug_predictions_file)) {
        drug_sensitivity <- read.csv(drug_predictions_file, row.names = 1, check.names = FALSE)
        # Save to output directory
        write.table(drug_sensitivity,
                    file.path(sensitivity_dir, "drug_sensitivity_scores.tsv"),
                    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
        cat(sprintf("  Predicted sensitivity for %d drugs x %d samples\n",
                    ncol(drug_sensitivity), nrow(drug_sensitivity)))
    } else {
        # Check if output went to a different location
        drug_sensitivity <- NULL
        cat("  WARNING: Drug predictions file not found at expected location\n")
    }

}, error = function(e) {
    cat(sprintf("  ERROR: oncoPredict failed: %s\n", conditionMessage(e)))
    drug_sensitivity <- NULL
})

# ============================================================================
# 2. Query DGIdb for Druggable Targets
# ============================================================================
cat("\n--- 2. DGIdb Druggable Target Query ---\n")

query_dgidb <- function(gene_list, batch_size = 100) {
    # DGIdb API endpoint
    base_url <- "https://dgidb.org/api/v2/interactions.json"
    all_interactions <- list()

    # Query in batches
    n_batches <- ceiling(length(gene_list) / batch_size)

    for (i in seq_len(n_batches)) {
        start_idx <- (i - 1) * batch_size + 1
        end_idx <- min(i * batch_size, length(gene_list))
        batch_genes <- gene_list[start_idx:end_idx]

        cat(sprintf("  Querying DGIdb batch %d/%d (%d genes)...\n", i, n_batches, length(batch_genes)))

        tryCatch({
            response <- GET(base_url,
                            query = list(genes = paste(batch_genes, collapse = ",")),
                            timeout(60))

            if (status_code(response) == 200) {
                result <- fromJSON(content(response, "text", encoding = "UTF-8"))

                if (!is.null(result$matchedTerms) && length(result$matchedTerms) > 0) {
                    for (term in seq_len(nrow(result$matchedTerms))) {
                        gene_name <- result$matchedTerms$searchTerm[term]
                        interactions <- result$matchedTerms$interactions[[term]]

                        if (!is.null(interactions) && nrow(interactions) > 0) {
                            interactions$gene <- gene_name
                            all_interactions[[length(all_interactions) + 1]] <- interactions
                        }
                    }
                }
            } else {
                cat(sprintf("    WARNING: DGIdb returned status %d for batch %d\n",
                            status_code(response), i))
            }

            # Rate limiting
            Sys.sleep(0.5)
        }, error = function(e) {
            cat(sprintf("    WARNING: DGIdb query failed for batch %d: %s\n",
                        i, conditionMessage(e)))
        })
    }

    if (length(all_interactions) > 0) {
        return(bind_rows(all_interactions))
    } else {
        return(data.frame())
    }
}

# Query DE genes against DGIdb
if (length(all_de_genes) > 0) {
    dgidb_results <- query_dgidb(all_de_genes)

    if (nrow(dgidb_results) > 0) {
        # Format results
        druggable_targets <- dgidb_results %>%
            select(gene, drugName, interactionType, score, drugChemblId, pmids) %>%
            distinct() %>%
            arrange(gene, desc(score))

        write.table(druggable_targets,
                    file.path(opt$output_dir, "druggable_targets", "druggable_targets.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        cat(sprintf("  Found %d drug-gene interactions for %d unique genes\n",
                    nrow(druggable_targets), length(unique(druggable_targets$gene))))

        # Summary of interaction types
        interaction_summary <- druggable_targets %>%
            count(interactionType, sort = TRUE)
        write.table(interaction_summary,
                    file.path(opt$output_dir, "druggable_targets", "interaction_type_summary.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
        cat("  No druggable targets found in DGIdb\n")
        druggable_targets <- data.frame()
    }
} else {
    cat("  No DE genes available for DGIdb query\n")
    druggable_targets <- data.frame()
}

# ============================================================================
# 3. Group Comparisons of Drug Sensitivity
# ============================================================================
cat("\n--- 3. Group Comparisons ---\n")

if (exists("drug_sensitivity") && !is.null(drug_sensitivity)) {

    comparison_results <- list()

    # Function to run Wilcoxon test between groups for all drugs
    compare_groups <- function(drug_mat, group_var, group_name) {
        groups <- unique(group_var)
        if (length(groups) != 2) return(NULL)

        results <- data.frame(
            drug = colnames(drug_mat),
            group1 = groups[1],
            group2 = groups[2],
            median_group1 = NA_real_,
            median_group2 = NA_real_,
            wilcox_pvalue = NA_real_,
            stringsAsFactors = FALSE
        )

        for (i in seq_len(ncol(drug_mat))) {
            g1_vals <- drug_mat[group_var == groups[1], i]
            g2_vals <- drug_mat[group_var == groups[2], i]

            results$median_group1[i] <- median(g1_vals, na.rm = TRUE)
            results$median_group2[i] <- median(g2_vals, na.rm = TRUE)

            if (length(g1_vals) >= 3 && length(g2_vals) >= 3) {
                test_result <- wilcox.test(g1_vals, g2_vals, exact = FALSE)
                results$wilcox_pvalue[i] <- test_result$p.value
            }
        }

        results$padj <- p.adjust(results$wilcox_pvalue, method = "BH")
        results$comparison <- group_name
        return(results)
    }

    # Compare by ancestry groups (if available)
    if (!is.null(ancestry)) {
        cat("  Comparing drug sensitivity by ancestry groups...\n")
        # Assign majority ancestry
        ancestry_aligned <- ancestry[rownames(drug_sensitivity), , drop = FALSE]
        ancestry_aligned <- ancestry_aligned[complete.cases(ancestry_aligned), , drop = FALSE]

        if (nrow(ancestry_aligned) > 0) {
            majority_ancestry <- colnames(ancestry_aligned)[apply(ancestry_aligned, 1, which.max)]
            names(majority_ancestry) <- rownames(ancestry_aligned)

            # For each pair of ancestry groups with sufficient samples
            anc_groups <- table(majority_ancestry)
            valid_groups <- names(anc_groups[anc_groups >= 5])

            if (length(valid_groups) >= 2) {
                for (i in 1:(length(valid_groups) - 1)) {
                    for (j in (i + 1):length(valid_groups)) {
                        g1 <- valid_groups[i]
                        g2 <- valid_groups[j]
                        mask <- majority_ancestry %in% c(g1, g2)
                        sub_drug <- drug_sensitivity[names(majority_ancestry)[mask], , drop = FALSE]
                        sub_groups <- majority_ancestry[mask]

                        comp_name <- paste0("ancestry_", g1, "_vs_", g2)
                        result <- compare_groups(sub_drug, sub_groups, comp_name)
                        if (!is.null(result)) {
                            comparison_results[[comp_name]] <- result
                        }
                    }
                }
            }
        }
    }

    # Compare by clinical groups from metadata
    clinical_cols <- c("condition", "group", "diagnosis", "risk_group", "subtype")
    for (clin_col in clinical_cols) {
        if (clin_col %in% colnames(metadata)) {
            cat(sprintf("  Comparing drug sensitivity by %s...\n", clin_col))
            meta_aligned <- metadata[rownames(drug_sensitivity), , drop = FALSE]
            meta_aligned <- meta_aligned[!is.na(meta_aligned[[clin_col]]), , drop = FALSE]

            grp_var <- meta_aligned[[clin_col]]
            grp_table <- table(grp_var)
            valid_grps <- names(grp_table[grp_table >= 5])

            if (length(valid_grps) == 2) {
                mask <- grp_var %in% valid_grps
                sub_drug <- drug_sensitivity[rownames(meta_aligned)[mask], , drop = FALSE]
                sub_groups <- grp_var[mask]
                comp_name <- paste0(clin_col, "_", valid_grps[1], "_vs_", valid_grps[2])
                result <- compare_groups(sub_drug, sub_groups, comp_name)
                if (!is.null(result)) {
                    comparison_results[[comp_name]] <- result
                }
            }
        }
    }

    # Save comparison results
    if (length(comparison_results) > 0) {
        all_comparisons <- bind_rows(comparison_results)
        write.table(all_comparisons,
                    file.path(opt$output_dir, "group_comparisons", "drug_sensitivity_comparisons.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        cat(sprintf("  Saved %d drug comparisons across %d group contrasts\n",
                    nrow(all_comparisons), length(comparison_results)))

        # Significant drugs per comparison
        sig_drugs <- all_comparisons %>%
            filter(padj < 0.05) %>%
            group_by(comparison) %>%
            summarise(n_significant_drugs = n(), .groups = "drop")
        cat("  Significant drugs per comparison:\n")
        print(sig_drugs)
    }
}

# ============================================================================
# Generate Plots
# ============================================================================
cat("\n--- Generating Plots ---\n")
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --- Drug sensitivity heatmap (top variable drugs) ---
if (exists("drug_sensitivity") && !is.null(drug_sensitivity)) {
    cat("  Creating drug sensitivity plots...\n")
    tryCatch({
        # Select top variable drugs
        drug_vars <- apply(drug_sensitivity, 2, var, na.rm = TRUE)
        top_drugs <- names(sort(drug_vars, decreasing = TRUE))[1:min(30, ncol(drug_sensitivity))]
        top_drug_mat <- as.matrix(drug_sensitivity[, top_drugs])

        # Boxplot of significant drugs by group
        if (length(comparison_results) > 0) {
            sig_drugs_list <- all_comparisons %>%
                filter(padj < 0.05) %>%
                arrange(padj) %>%
                slice_head(n = 10) %>%
                pull(drug)

            if (length(sig_drugs_list) > 0) {
                plot_drugs <- sig_drugs_list[sig_drugs_list %in% colnames(drug_sensitivity)]
                if (length(plot_drugs) > 0) {
                    # Get first comparison group info
                    first_comp <- comparison_results[[1]]
                    comp_name <- unique(first_comp$comparison)[1]

                    drug_long <- drug_sensitivity[, plot_drugs, drop = FALSE] %>%
                        rownames_to_column("sample") %>%
                        pivot_longer(cols = -sample, names_to = "drug", values_to = "sensitivity")

                    # Add group information
                    if (!is.null(ancestry)) {
                        ancestry_aligned <- ancestry[drug_long$sample, , drop = FALSE]
                        if (nrow(ancestry_aligned) > 0) {
                            drug_long$group <- colnames(ancestry)[apply(
                                ancestry[drug_long$sample, , drop = FALSE], 1, which.max
                            )]
                        }
                    }

                    if ("group" %in% colnames(drug_long)) {
                        p_box <- ggplot(drug_long, aes(x = drug, y = sensitivity, fill = group)) +
                            geom_boxplot(outlier.size = 0.5) +
                            theme_minimal() +
                            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
                            labs(
                                title = "Drug Sensitivity Scores by Group (Top Significant)",
                                x = "Drug",
                                y = "Predicted IC50 (log)",
                                fill = "Group"
                            )
                        ggsave(file.path(plot_dir, "drug_sensitivity_boxplot.pdf"),
                               p_box, width = 12, height = 6)
                    }
                }
            }
        }

        # Volcano-like plot of drug sensitivity differences
        if (length(comparison_results) > 0) {
            for (comp_name in names(comparison_results)) {
                comp_df <- comparison_results[[comp_name]]
                comp_df$neg_log10_p <- -log10(comp_df$wilcox_pvalue)
                comp_df$diff_median <- comp_df$median_group1 - comp_df$median_group2
                comp_df$significant <- comp_df$padj < 0.05

                p_volcano <- ggplot(comp_df, aes(x = diff_median, y = neg_log10_p,
                                                  color = significant)) +
                    geom_point(alpha = 0.6) +
                    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
                    theme_minimal() +
                    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
                    labs(
                        title = paste("Drug Sensitivity Differences:", comp_name),
                        x = "Difference in Median IC50",
                        y = "-log10(p-value)",
                        color = "Significant (FDR<0.05)"
                    )
                ggsave(file.path(plot_dir, paste0("volcano_", comp_name, ".pdf")),
                       p_volcano, width = 10, height = 7)
            }
        }
        cat("    Drug sensitivity plots saved.\n")
    }, error = function(e) {
        cat(sprintf("    WARNING: Drug sensitivity plots failed: %s\n", conditionMessage(e)))
    })
}

# --- Druggable targets barplot ---
if (nrow(druggable_targets) > 0) {
    cat("  Creating druggable targets plot...\n")
    tryCatch({
        # Top genes by number of drug interactions
        gene_drug_counts <- druggable_targets %>%
            count(gene, sort = TRUE) %>%
            slice_head(n = 30)

        p_targets <- ggplot(gene_drug_counts, aes(x = reorder(gene, n), y = n)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            coord_flip() +
            theme_minimal() +
            labs(
                title = "Top Druggable DE Genes (by Number of Drug Interactions)",
                x = "Gene",
                y = "Number of Drug Interactions"
            )
        ggsave(file.path(plot_dir, "druggable_targets_barplot.pdf"),
               p_targets, width = 10, height = 8)
        cat("    Druggable targets plot saved.\n")
    }, error = function(e) {
        cat(sprintf("    WARNING: Druggable targets plot failed: %s\n", conditionMessage(e)))
    })
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Pharmacogenomics Analysis Complete ===\n")
cat(sprintf("Output directory: %s\n", opt$output_dir))
cat("Key outputs:\n")
cat(sprintf("  - Drug sensitivity scores: %s\n",
            file.path(opt$output_dir, "drug_sensitivity", "drug_sensitivity_scores.tsv")))
cat(sprintf("  - Druggable targets: %s\n",
            file.path(opt$output_dir, "druggable_targets", "druggable_targets.tsv")))
cat(sprintf("  - Group comparisons: %s\n",
            file.path(opt$output_dir, "group_comparisons", "drug_sensitivity_comparisons.tsv")))
cat(sprintf("  - Plots: %s\n", plot_dir))
