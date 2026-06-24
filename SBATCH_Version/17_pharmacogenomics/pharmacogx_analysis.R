#!/usr/bin/env Rscript
# ============================================================================
# 17 - PharmacoGx Drug Sensitivity Analysis
# ============================================================================
# Downloads PharmacoSets (GDSC2, CTRPv2, gCSI), trains ridge regression models
# for drug sensitivity prediction using cell line expression + drug response
# data, predicts per-sample drug sensitivity scores from tumor RNA-seq, and
# compares predictions across ancestry groups (Kruskal-Wallis).
#
# Required packages: PharmacoGx, glmnet, ComplexHeatmap, ggplot2, tidyverse,
#                    optparse
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(PharmacoGx)
    library(glmnet)
    library(tidyverse)
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
    make_option("--normalized_counts", type = "character",
                help = "Path to normalized count matrix (genes x samples TSV)"),
    make_option("--metadata", type = "character",
                help = "Path to sample metadata CSV"),
    make_option("--ancestry_proportions", type = "character", default = "none",
                help = "Path to ancestry proportions TSV [default: none]"),
    make_option("--output_dir", type = "character",
                help = "Output directory"),
    make_option("--datasets", type = "character", default = "GDSC2,CTRPv2,gCSI",
                help = "Comma-separated PharmacoSet names [default: %default]"),
    make_option("--threads", type = "integer", default = 8,
                help = "Number of threads [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$normalized_counts) || is.null(opt$metadata) ||
    is.null(opt$output_dir)) {
    stop("Required: --normalized_counts, --metadata, --output_dir")
}

cat("=== PharmacoGx Drug Sensitivity Analysis ===\n")
cat("Normalized counts:", opt$normalized_counts, "\n")
cat("Metadata:", opt$metadata, "\n")
cat("Ancestry proportions:", opt$ancestry_proportions, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Datasets:", opt$datasets, "\n")
cat("Threads:", opt$threads, "\n\n")

# ============================================================================
# Create output directories
# ============================================================================
pred_dir <- file.path(opt$output_dir, "predictions")
comp_dir <- file.path(opt$output_dir, "comparisons")
model_dir <- file.path(opt$output_dir, "models")
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Load data
# ============================================================================
cat("Loading data...\n")

counts <- read.table(opt$normalized_counts, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE)
metadata <- read.csv(opt$metadata, header = TRUE, stringsAsFactors = FALSE)

cat(sprintf("  Count matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  Metadata: %d samples\n", nrow(metadata)))

# Align metadata with counts
sample_ids <- colnames(counts)
metadata <- metadata %>% filter(metadata[[1]] %in% sample_ids)
rownames(metadata) <- metadata[[1]]

# Load ancestry proportions
ancestry <- NULL
majority_ancestry <- NULL
if (opt$ancestry_proportions != "none" && file.exists(opt$ancestry_proportions)) {
    ancestry <- read.table(opt$ancestry_proportions, header = TRUE, sep = "\t",
                           row.names = 1, check.names = FALSE)
    # Align with count matrix samples
    common_samples <- intersect(rownames(ancestry), sample_ids)
    if (length(common_samples) > 0) {
        ancestry <- ancestry[common_samples, , drop = FALSE]
        majority_ancestry <- colnames(ancestry)[apply(ancestry, 1, which.max)]
        names(majority_ancestry) <- common_samples
        cat(sprintf("  Ancestry: %d samples, %d groups (%s)\n",
                    length(common_samples), length(unique(majority_ancestry)),
                    paste(unique(majority_ancestry), collapse = ", ")))
    }
}

# ============================================================================
# Parse dataset list
# ============================================================================
dataset_names <- trimws(unlist(strsplit(opt$datasets, ",")))
cat(sprintf("\nWill process %d PharmacoSets: %s\n",
            length(dataset_names), paste(dataset_names, collapse = ", ")))

# ============================================================================
# Process each PharmacoSet
# ============================================================================
all_predictions <- list()

for (ds_name in dataset_names) {
    cat(sprintf("\n--- Processing PharmacoSet: %s ---\n", ds_name))

    tryCatch({
        # Download PharmacoSet
        cat(sprintf("  Downloading %s PharmacoSet...\n", ds_name))
        pset <- downloadPSet(ds_name, saveDir = file.path(opt$output_dir, "psets"))

        # Extract expression data from cell lines
        cat("  Extracting cell line expression data...\n")
        cl_expr <- summarizeMolecularProfiles(pset, mDataType = "rna",
                                               cell.lines = cellNames(pset),
                                               features = featureNames(pset, "rna"))
        cl_expr_mat <- assay(cl_expr)

        # Extract drug sensitivity (AAC = area above the curve)
        cat("  Extracting drug sensitivity data...\n")
        drug_resp <- summarizeSensitivityProfiles(pset, sensitivity.measure = "aac_recomputed",
                                                   summary.stat = "median")

        cat(sprintf("  Cell line expression: %d genes x %d cell lines\n",
                    nrow(cl_expr_mat), ncol(cl_expr_mat)))
        cat(sprintf("  Drug response: %d cell lines x %d drugs\n",
                    nrow(drug_resp), ncol(drug_resp)))

        # Find common genes between cell lines and tumor samples
        common_genes <- intersect(rownames(cl_expr_mat), rownames(counts))
        cat(sprintf("  Common genes: %d\n", length(common_genes)))

        if (length(common_genes) < 100) {
            cat(sprintf("  WARNING: Too few common genes for %s, skipping\n", ds_name))
            next
        }

        # Align expression matrices to common genes
        cl_train <- t(cl_expr_mat[common_genes, , drop = FALSE])
        tumor_test <- t(counts[common_genes, , drop = FALSE])

        # Find common cell lines between expression and drug response
        common_cls <- intersect(rownames(cl_train), rownames(drug_resp))
        cl_train <- cl_train[common_cls, , drop = FALSE]
        drug_resp_aligned <- drug_resp[common_cls, , drop = FALSE]

        cat(sprintf("  Aligned: %d cell lines, %d genes, %d drugs\n",
                    nrow(cl_train), ncol(cl_train), ncol(drug_resp_aligned)))

        # -----------------------------------------------------------------
        # Train ridge regression models for each drug
        # -----------------------------------------------------------------
        cat("  Training ridge regression models...\n")
        predictions <- matrix(NA, nrow = nrow(tumor_test), ncol = ncol(drug_resp_aligned))
        rownames(predictions) <- rownames(tumor_test)
        colnames(predictions) <- colnames(drug_resp_aligned)
        model_count <- 0

        for (drug_idx in seq_len(ncol(drug_resp_aligned))) {
            drug_name <- colnames(drug_resp_aligned)[drug_idx]
            y <- drug_resp_aligned[, drug_idx]
            valid <- !is.na(y)

            if (sum(valid) < 20) next

            tryCatch({
                # Fit ridge regression (alpha = 0) with cross-validation
                cv_fit <- cv.glmnet(
                    x = cl_train[valid, , drop = FALSE],
                    y = y[valid],
                    alpha = 0,
                    nfolds = 5,
                    type.measure = "mse"
                )

                # Predict on tumor samples
                pred <- predict(cv_fit, newx = tumor_test, s = "lambda.min")
                predictions[, drug_idx] <- pred[, 1]
                model_count <- model_count + 1
            }, error = function(e) NULL)
        }

        # Remove drugs with no predictions
        valid_drugs <- colSums(!is.na(predictions)) > 0
        predictions <- predictions[, valid_drugs, drop = FALSE]

        cat(sprintf("  Successfully trained models for %d / %d drugs\n",
                    model_count, ncol(drug_resp_aligned)))

        # Save predictions
        pred_file <- file.path(pred_dir, paste0(ds_name, "_drug_sensitivity.tsv"))
        write.table(predictions, pred_file,
                    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
        cat(sprintf("  Predictions saved: %s\n", pred_file))

        all_predictions[[ds_name]] <- predictions

        # -----------------------------------------------------------------
        # Ancestry group comparisons (Kruskal-Wallis)
        # -----------------------------------------------------------------
        if (!is.null(majority_ancestry) && ncol(predictions) > 0) {
            cat("  Comparing drug sensitivity across ancestry groups...\n")

            # Filter to samples with ancestry and predictions
            common_pred_anc <- intersect(rownames(predictions), names(majority_ancestry))

            if (length(common_pred_anc) >= 10) {
                pred_sub <- predictions[common_pred_anc, , drop = FALSE]
                anc_sub <- majority_ancestry[common_pred_anc]

                # Only test if at least 2 groups with >= 3 samples
                anc_table <- table(anc_sub)
                valid_groups <- names(anc_table[anc_table >= 3])

                if (length(valid_groups) >= 2) {
                    mask <- anc_sub %in% valid_groups
                    pred_sub <- pred_sub[mask, , drop = FALSE]
                    anc_sub <- anc_sub[mask]

                    kw_results <- data.frame(
                        drug = colnames(pred_sub),
                        kruskal_wallis_pvalue = NA_real_,
                        stringsAsFactors = FALSE
                    )

                    # Add median per group
                    for (grp in valid_groups) {
                        kw_results[[paste0("median_", grp)]] <- NA_real_
                    }

                    for (i in seq_len(ncol(pred_sub))) {
                        vals <- pred_sub[, i]
                        if (sum(!is.na(vals)) >= length(valid_groups) * 3) {
                            test_res <- kruskal.test(vals ~ factor(anc_sub))
                            kw_results$kruskal_wallis_pvalue[i] <- test_res$p.value

                            for (grp in valid_groups) {
                                kw_results[[paste0("median_", grp)]][i] <-
                                    median(vals[anc_sub == grp], na.rm = TRUE)
                            }
                        }
                    }

                    kw_results$padj <- p.adjust(kw_results$kruskal_wallis_pvalue,
                                                method = "BH")
                    kw_results <- kw_results %>% arrange(padj)

                    comp_file <- file.path(comp_dir,
                                           paste0(ds_name, "_ancestry_comparison.tsv"))
                    write.table(kw_results, comp_file,
                                sep = "\t", quote = FALSE, row.names = FALSE)

                    n_sig <- sum(kw_results$padj < 0.05, na.rm = TRUE)
                    cat(sprintf("  %d drugs with significant ancestry differences (FDR < 0.05)\n",
                                n_sig))
                }
            }
        }

        # -----------------------------------------------------------------
        # Drug sensitivity heatmap
        # -----------------------------------------------------------------
        cat("  Generating drug sensitivity heatmap...\n")
        tryCatch({
            # Select top variable drugs for heatmap
            drug_vars <- apply(predictions, 2, var, na.rm = TRUE)
            top_drugs <- names(sort(drug_vars, decreasing = TRUE))[
                1:min(50, ncol(predictions))
            ]
            hm_mat <- t(scale(predictions[, top_drugs, drop = FALSE]))

            # Build annotation
            ha_list <- list()
            if (!is.null(majority_ancestry)) {
                common_hm <- intersect(colnames(hm_mat), names(majority_ancestry))
                if (length(common_hm) > 0) {
                    hm_mat <- hm_mat[, common_hm, drop = FALSE]
                    ha_list[["Ancestry"]] <- majority_ancestry[common_hm]
                }
            }

            col_ha <- NULL
            if (length(ha_list) > 0) {
                col_ha <- HeatmapAnnotation(
                    Ancestry = ha_list[["Ancestry"]],
                    annotation_name_side = "left"
                )
            }

            col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

            pdf(file.path(plot_dir, paste0(ds_name, "_drug_sensitivity_heatmap.pdf")),
                width = 14, height = 10)
            ht <- Heatmap(
                hm_mat,
                name = "Scaled\nSensitivity",
                col = col_fun,
                top_annotation = col_ha,
                show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 7),
                column_title = paste("Drug Sensitivity Predictions -", ds_name),
                clustering_distance_rows = "euclidean",
                clustering_distance_columns = "euclidean",
                clustering_method_rows = "ward.D2",
                clustering_method_columns = "ward.D2"
            )
            draw(ht)
            dev.off()

            cat("    Heatmap saved.\n")
        }, error = function(e) {
            cat(sprintf("    WARNING: Heatmap generation failed: %s\n",
                        conditionMessage(e)))
        })

    }, error = function(e) {
        cat(sprintf("  ERROR processing %s: %s\n", ds_name, conditionMessage(e)))
    })
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== PharmacoGx Analysis Complete ===\n")
cat(sprintf("Output directory: %s\n", opt$output_dir))
cat("Key outputs:\n")
cat(sprintf("  - Drug sensitivity predictions: %s\n", pred_dir))
cat(sprintf("  - Ancestry comparisons: %s\n", comp_dir))
cat(sprintf("  - Plots: %s\n", plot_dir))
for (ds in names(all_predictions)) {
    cat(sprintf("  - %s: %d drugs x %d samples\n",
                ds, ncol(all_predictions[[ds]]), nrow(all_predictions[[ds]])))
}
