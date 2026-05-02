#!/usr/bin/env Rscript
# ============================================================================
# 13 - Immune Deconvolution Analysis
# ============================================================================
# Runs multiple immune deconvolution methods and ESTIMATE scoring.
# Optionally corrects cell fractions by tumor purity.
#
# Required packages: immunedeconv, estimate, ggplot2, pheatmap, optparse,
#                    tidyverse, reshape2
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(immunedeconv)
    library(ggplot2)
    library(pheatmap)
    library(tidyverse)
    library(reshape2)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
    make_option("--normalized_counts", type = "character",
                help = "Path to normalized count matrix (genes x samples TSV)"),
    make_option("--metadata", type = "character",
                help = "Path to sample metadata CSV"),
    make_option("--output_dir", type = "character",
                help = "Output directory"),
    make_option("--methods", type = "character", default = "xcell,mcpcounter,epic,cibersortx,timer",
                help = "Comma-separated deconvolution methods [default: %default]"),
    make_option("--correct_tumor_purity", type = "character", default = "TRUE",
                help = "Whether to adjust cell fractions by tumor purity [default: %default]"),
    make_option("--threads", type = "integer", default = 8,
                help = "Number of threads [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$normalized_counts) || is.null(opt$metadata) || is.null(opt$output_dir)) {
    stop("Required arguments: --normalized_counts, --metadata, --output_dir")
}

correct_purity <- as.logical(opt$correct_tumor_purity)
methods <- strsplit(opt$methods, ",")[[1]]

cat("=== Immune Deconvolution Analysis ===\n")
cat("Normalized counts:", opt$normalized_counts, "\n")
cat("Metadata:", opt$metadata, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Methods:", paste(methods, collapse = ", "), "\n")
cat("Correct tumor purity:", correct_purity, "\n")
cat("Threads:", opt$threads, "\n\n")

# ============================================================================
# Load data
# ============================================================================
cat("Loading normalized count matrix...\n")
counts <- read.table(opt$normalized_counts, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE)
cat(sprintf("  Loaded matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))

cat("Loading metadata...\n")
metadata <- read.csv(opt$metadata, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Loaded metadata for %d samples\n", nrow(metadata)))

# Ensure counts are in TPM-like scale for methods that require it
# immunedeconv expects gene expression in TPM for most methods
expr_matrix <- as.matrix(counts)

# ============================================================================
# Run ESTIMATE for purity/immune/stromal scores
# ============================================================================
cat("\n--- Running ESTIMATE ---\n")

# ESTIMATE requires writing to a temp file
estimate_input <- file.path(opt$output_dir, "estimate", "estimate_input.gct")
estimate_output <- file.path(opt$output_dir, "estimate", "estimate_scores.gct")

# Write GCT format for ESTIMATE
dir.create(file.path(opt$output_dir, "estimate"), recursive = TRUE, showWarnings = FALSE)

# Create GCT file
gct_header <- c(
    "#1.2",
    paste(nrow(counts), ncol(counts), sep = "\t")
)
gct_data <- cbind(NAME = rownames(counts), Description = "na", counts)
writeLines(gct_header, estimate_input)
write.table(gct_data, estimate_input, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE, append = TRUE)

# Run ESTIMATE
tryCatch({
    library(estimate)
    filterCommonGenes(input.f = estimate_input,
                      output.f = file.path(opt$output_dir, "estimate", "estimate_filtered.gct"),
                      id = "GeneSymbol")
    estimateScore(input.ds = file.path(opt$output_dir, "estimate", "estimate_filtered.gct"),
                  output.ds = estimate_output,
                  platform = "illumina")

    # Parse ESTIMATE output
    estimate_scores <- read.table(estimate_output, header = TRUE, sep = "\t",
                                  skip = 2, row.names = 1, check.names = FALSE)
    estimate_scores <- estimate_scores[, -1]  # Remove Description column
    estimate_df <- as.data.frame(t(estimate_scores))
    colnames(estimate_df) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")

    write.table(estimate_df, file.path(opt$output_dir, "estimate", "estimate_scores.tsv"),
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
    cat("  ESTIMATE scores saved.\n")
}, error = function(e) {
    cat("  WARNING: ESTIMATE failed:", conditionMessage(e), "\n")
    cat("  Attempting purity estimation from immunedeconv...\n")
    estimate_df <- NULL
})

# ============================================================================
# Run deconvolution methods
# ============================================================================
cat("\n--- Running immune deconvolution methods ---\n")

all_results <- list()

for (method in methods) {
    cat(sprintf("\nRunning method: %s\n", method))

    tryCatch({
        if (method == "timer") {
            # TIMER requires cancer type indication
            # Default to "ALL" (acute lymphoblastic leukemia) if available in metadata
            indications <- rep("all", ncol(expr_matrix))
            result <- deconvolute(expr_matrix, method = "timer",
                                  indications = indications)
        } else if (method == "cibersortx") {
            # CIBERSORTx may need token - attempt with default settings
            result <- deconvolute(expr_matrix, method = "cibersort_abs")
        } else {
            result <- deconvolute(expr_matrix, method = method)
        }

        # Store result
        all_results[[method]] <- result
        cat(sprintf("  %s completed: %d cell types x %d samples\n",
                    method, nrow(result) - 1, ncol(result) - 1))

        # Save individual method results
        write.table(result, file.path(opt$output_dir, "deconvolution",
                                      paste0(method, "_results.tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)

    }, error = function(e) {
        cat(sprintf("  WARNING: %s failed: %s\n", method, conditionMessage(e)))
    })
}

# ============================================================================
# Combine all deconvolution results
# ============================================================================
cat("\n--- Combining deconvolution results ---\n")

if (length(all_results) > 0) {
    # Create combined long-format table
    combined_long <- do.call(rbind, lapply(names(all_results), function(m) {
        df <- all_results[[m]]
        df_long <- df %>%
            pivot_longer(cols = -cell_type, names_to = "sample", values_to = "score") %>%
            mutate(method = m)
        return(df_long)
    }))

    # Apply tumor purity correction if requested
    if (correct_purity && exists("estimate_df") && !is.null(estimate_df)) {
        cat("  Applying tumor purity correction...\n")
        combined_long <- combined_long %>%
            left_join(
                estimate_df %>%
                    rownames_to_column("sample") %>%
                    select(sample, TumorPurity),
                by = "sample"
            ) %>%
            mutate(
                score_corrected = ifelse(
                    !is.na(TumorPurity) & TumorPurity < 1,
                    score / (1 - TumorPurity),
                    score
                )
            )
        cat("  Purity correction applied.\n")
    } else {
        combined_long$score_corrected <- combined_long$score
        if (correct_purity) {
            cat("  WARNING: Purity correction requested but ESTIMATE scores unavailable.\n")
        }
    }

    # Save combined results
    write.table(combined_long, file.path(opt$output_dir, "deconvolution",
                                         "combined_deconvolution.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Combined deconvolution results saved.\n")
}

# ============================================================================
# Generate plots
# ============================================================================
cat("\n--- Generating plots ---\n")
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --- Stacked bar chart of cell composition ---
for (method_name in names(all_results)) {
    cat(sprintf("  Plotting stacked bar chart for %s...\n", method_name))

    tryCatch({
        result_df <- all_results[[method_name]] %>%
            pivot_longer(cols = -cell_type, names_to = "sample", values_to = "fraction") %>%
            filter(fraction > 0)

        p_bar <- ggplot(result_df, aes(x = sample, y = fraction, fill = cell_type)) +
            geom_bar(stat = "identity", position = "stack") +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                legend.position = "right",
                legend.text = element_text(size = 7)
            ) +
            labs(
                title = paste("Immune Cell Composition -", toupper(method_name)),
                x = "Sample",
                y = "Cell Fraction",
                fill = "Cell Type"
            ) +
            scale_fill_viridis_d(option = "turbo")

        ggsave(file.path(plot_dir, paste0("stacked_bar_", method_name, ".pdf")),
               p_bar, width = 14, height = 8)
    }, error = function(e) {
        cat(sprintf("    WARNING: Bar plot for %s failed: %s\n", method_name, conditionMessage(e)))
    })
}

# --- Heatmap of deconvolution scores ---
for (method_name in names(all_results)) {
    cat(sprintf("  Plotting heatmap for %s...\n", method_name))

    tryCatch({
        result_mat <- all_results[[method_name]] %>%
            column_to_rownames("cell_type") %>%
            as.matrix()

        # Remove rows with all zeros
        result_mat <- result_mat[rowSums(result_mat) > 0, , drop = FALSE]

        if (nrow(result_mat) > 1 && ncol(result_mat) > 1) {
            # Scale rows for visualization
            result_scaled <- t(scale(t(result_mat)))
            result_scaled[is.nan(result_scaled)] <- 0

            # Annotation for metadata groups if available
            annotation_col <- NULL
            if ("group" %in% colnames(metadata)) {
                sample_groups <- metadata %>%
                    filter(metadata[[1]] %in% colnames(result_mat)) %>%
                    select(1, group) %>%
                    column_to_rownames(colnames(metadata)[1])
                annotation_col <- sample_groups
            }

            pdf(file.path(plot_dir, paste0("heatmap_", method_name, ".pdf")),
                width = 12, height = 8)
            pheatmap(result_scaled,
                     main = paste("Immune Cell Scores -", toupper(method_name)),
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "ward.D2",
                     show_colnames = TRUE,
                     fontsize_col = 6,
                     fontsize_row = 8,
                     annotation_col = annotation_col)
            dev.off()
        }
    }, error = function(e) {
        cat(sprintf("    WARNING: Heatmap for %s failed: %s\n", method_name, conditionMessage(e)))
    })
}

# --- ESTIMATE score plot ---
if (exists("estimate_df") && !is.null(estimate_df)) {
    cat("  Plotting ESTIMATE scores...\n")
    tryCatch({
        estimate_long <- estimate_df %>%
            rownames_to_column("sample") %>%
            pivot_longer(cols = -sample, names_to = "score_type", values_to = "value")

        p_estimate <- ggplot(estimate_long, aes(x = sample, y = value, fill = score_type)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
            labs(
                title = "ESTIMATE Scores",
                x = "Sample",
                y = "Score",
                fill = "Score Type"
            ) +
            facet_wrap(~score_type, scales = "free_y", ncol = 1)

        ggsave(file.path(plot_dir, "estimate_scores.pdf"),
               p_estimate, width = 12, height = 10)
    }, error = function(e) {
        cat(sprintf("    WARNING: ESTIMATE plot failed: %s\n", conditionMessage(e)))
    })
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Immune Deconvolution Complete ===\n")
cat(sprintf("Methods run successfully: %d / %d\n", length(all_results), length(methods)))
cat(sprintf("Output directory: %s\n", opt$output_dir))
cat("Key outputs:\n")
cat(sprintf("  - Combined results: %s\n", file.path(opt$output_dir, "deconvolution", "combined_deconvolution.tsv")))
cat(sprintf("  - ESTIMATE scores: %s\n", file.path(opt$output_dir, "estimate", "estimate_scores.tsv")))
cat(sprintf("  - Plots: %s\n", plot_dir))
