#!/usr/bin/env Rscript

# ============================================================================
# Expression Outlier Detection
# Identifies potential outlier samples using three complementary approaches:
#   1. PCA-based: samples distant from centroid in PC space
#   2. Correlation-based: samples with low pairwise correlation to others
#   3. Library size: samples with extreme total counts
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(stats)
})

# ---- Parse arguments ----
option_list <- list(
    make_option("--counts", type = "character",
                help = "Count matrix TSV (genes x samples, raw or normalized)"),
    make_option("--metadata", type = "character",
                help = "Sample metadata TSV with 'sample_id' column"),
    make_option("--output-dir", type = "character", default = "outlier_results",
                help = "Output directory [default: %default]"),
    make_option("--sd-threshold", type = "double", default = 3,
                help = "Standard deviations from PCA centroid to flag outlier [default: %default]"),
    make_option("--min-correlation", type = "double", default = 0.8,
                help = "Minimum median pairwise correlation threshold [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate inputs ----
if (is.null(opt$counts)) {
    stop("--counts argument is required.")
}

if (!file.exists(opt$counts)) {
    stop(paste("Counts file not found:", opt$counts))
}

if (!is.null(opt$metadata) && !file.exists(opt$metadata)) {
    stop(paste("Metadata file not found:", opt$metadata))
}

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

sd_threshold <- opt$`sd-threshold`
min_correlation <- opt$`min-correlation`

cat("Parameters:\n")
cat("  SD threshold:", sd_threshold, "\n")
cat("  Min correlation:", min_correlation, "\n")

# ---- Load data ----
cat("Loading count matrix:", opt$counts, "\n")
counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)

metadata <- NULL
if (!is.null(opt$metadata)) {
    cat("Loading metadata:", opt$metadata, "\n")
    metadata <- read.delim(opt$metadata, check.names = FALSE)
    rownames(metadata) <- metadata$sample_id

    common_samples <- intersect(colnames(counts), rownames(metadata))
    if (length(common_samples) > 0) {
        counts <- counts[, common_samples, drop = FALSE]
        metadata <- metadata[common_samples, , drop = FALSE]
    }
}

n_samples <- ncol(counts)
n_genes <- nrow(counts)
cat("Samples:", n_samples, "\n")
cat("Genes:", n_genes, "\n")

if (n_samples < 3) {
    findings <- c(
        "# Expression Outlier Detection",
        "",
        "## Status: SKIPPED",
        "",
        paste0("Too few samples (n=", n_samples, ") for outlier detection. ",
               "At least 3 samples are required."),
        ""
    )
    writeLines(findings, file.path(opt$`output-dir`, "expression_outlier_findings.md"))
    cat("Too few samples. Findings written.\n")
    quit(status = 0)
}

# ---- Prepare data: log2-transform and filter ----
# Filter low-expression genes
gene_means <- rowMeans(as.matrix(counts))
keep_genes <- gene_means > 1
counts_filtered <- counts[keep_genes, , drop = FALSE]
cat("Genes after filtering (mean > 1):", nrow(counts_filtered), "\n")

if (nrow(counts_filtered) < 100) {
    warning("Fewer than 100 genes passed filtering. Results may be unreliable.")
}

# Log2 transform with pseudocount
log2_counts <- log2(as.matrix(counts_filtered) + 1)

# ---- 1. PCA-based outlier detection ----
cat("\n==== PCA-based outlier detection ====\n")

# Center and scale for PCA
log2_scaled <- t(scale(t(log2_counts)))
# Remove genes with zero variance
zero_var <- apply(log2_scaled, 1, function(x) all(is.na(x)) || var(x, na.rm = TRUE) == 0)
log2_scaled <- log2_scaled[!zero_var, ]
log2_scaled[is.na(log2_scaled)] <- 0

pca_result <- prcomp(t(log2_scaled), center = TRUE, scale. = FALSE)
pca_scores <- pca_result$x

# Variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100

# Compute Euclidean distance from centroid in top PCs (capturing >= 50% variance)
cumvar <- cumsum(var_explained)
n_pcs <- max(2, min(which(cumvar >= 50), n_samples - 1))
n_pcs <- min(n_pcs, ncol(pca_scores))
cat("Using top", n_pcs, "PCs (cumulative variance:",
    round(cumvar[n_pcs], 1), "%)\n")

pc_subset <- pca_scores[, 1:n_pcs, drop = FALSE]
centroid <- colMeans(pc_subset)
distances <- sqrt(rowSums(sweep(pc_subset, 2, centroid)^2))

# Flag outliers
dist_mean <- mean(distances)
dist_sd <- sd(distances)
pca_outlier <- distances > (dist_mean + sd_threshold * dist_sd)

cat("PCA outliers (>", sd_threshold, "SD from centroid):",
    sum(pca_outlier), "\n")

# ---- 2. Correlation-based outlier detection ----
cat("\n==== Correlation-based outlier detection ====\n")

# Compute pairwise Spearman correlation on a subset of variable genes for efficiency
gene_vars <- apply(log2_counts, 1, var)
top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(5000, length(gene_vars))]
cor_matrix <- cor(log2_counts[top_var_genes, ], method = "spearman", use = "pairwise.complete.obs")

# Median pairwise correlation for each sample (excluding self)
median_cor <- sapply(1:n_samples, function(i) {
    cors <- cor_matrix[i, -i]
    median(cors, na.rm = TRUE)
})
names(median_cor) <- colnames(counts)

cor_outlier <- median_cor < min_correlation
cat("Correlation outliers (median r <", min_correlation, "):",
    sum(cor_outlier), "\n")

# ---- 3. Library size check ----
cat("\n==== Library size check ====\n")

lib_sizes <- colSums(as.matrix(counts))
log_lib <- log10(lib_sizes + 1)
lib_mean <- mean(log_lib)
lib_sd <- sd(log_lib)

# Flag samples with library size beyond sd_threshold SDs (on log scale)
lib_outlier_high <- log_lib > (lib_mean + sd_threshold * lib_sd)
lib_outlier_low <- log_lib < (lib_mean - sd_threshold * lib_sd)
lib_outlier <- lib_outlier_high | lib_outlier_low

cat("Library size outliers:", sum(lib_outlier), "\n")
cat("  High:", sum(lib_outlier_high), "| Low:", sum(lib_outlier_low), "\n")
cat("  Library size range:", sprintf("%.2e - %.2e", min(lib_sizes), max(lib_sizes)), "\n")
cat("  Median library size:", sprintf("%.2e", median(lib_sizes)), "\n")

# ---- Combine results ----
results_df <- data.frame(
    sample_id = colnames(counts),
    pca_distance = round(distances, 4),
    pca_distance_zscore = round((distances - dist_mean) / dist_sd, 4),
    pca_outlier = pca_outlier,
    median_pairwise_correlation = round(median_cor, 4),
    correlation_outlier = cor_outlier,
    library_size = lib_sizes,
    log10_library_size = round(log_lib, 4),
    library_size_zscore = round((log_lib - lib_mean) / lib_sd, 4),
    library_size_outlier = lib_outlier,
    n_flags = as.integer(pca_outlier) + as.integer(cor_outlier) + as.integer(lib_outlier),
    stringsAsFactors = FALSE
)

# Sort by number of flags (most flagged first)
results_df <- results_df[order(-results_df$n_flags, -results_df$pca_distance_zscore), ]

write.table(results_df, file.path(opt$`output-dir`, "outlier_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nResults written to outlier_results.tsv\n")

# ---- Generate PCA outlier plot ----
tryCatch({
    pdf(file.path(opt$`output-dir`, "pca_outlier_plot.pdf"), width = 9, height = 7)

    plot_df <- data.frame(
        PC1 = pca_scores[, 1],
        PC2 = pca_scores[, 2],
        sample_id = colnames(counts),
        outlier = ifelse(pca_outlier, "Outlier", "Normal"),
        n_flags = results_df[match(colnames(counts), results_df$sample_id), "n_flags"]
    )

    p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = outlier)) +
        geom_point(size = 3, alpha = 0.7) +
        scale_color_manual(values = c("Normal" = "steelblue", "Outlier" = "red")) +
        labs(
            title = "PCA-based Outlier Detection",
            subtitle = paste0("SD threshold: ", sd_threshold,
                             " | Outliers: ", sum(pca_outlier), "/", n_samples),
            x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
            y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
            color = "Status"
        ) +
        theme_bw(base_size = 12)

    # Label outlier samples
    if (sum(pca_outlier) > 0) {
        outlier_df <- plot_df[plot_df$outlier == "Outlier", ]
        p <- p + geom_text(data = outlier_df,
                            aes(label = sample_id),
                            vjust = -1, size = 2.5, color = "red")
    }

    print(p)
    dev.off()
    cat("PCA outlier plot written.\n")
}, error = function(e) {
    cat("PCA outlier plot failed:", e$message, "\n")
})

# ---- Generate correlation heatmap ----
tryCatch({
    pdf(file.path(opt$`output-dir`, "correlation_heatmap.pdf"),
        width = max(8, n_samples * 0.3), height = max(7, n_samples * 0.3))

    # Annotate outlier status
    sample_colors <- ifelse(cor_outlier, "red", "black")
    names(sample_colors) <- colnames(counts)

    # If too many samples, cluster and show dendrogram only
    if (n_samples <= 100) {
        hc <- hclust(as.dist(1 - cor_matrix), method = "ward.D2")

        # Use base R heatmap for robustness (no extra dependencies)
        heatmap(cor_matrix,
                Rowv = as.dendrogram(hc),
                Colv = as.dendrogram(hc),
                col = colorRampPalette(c("#313695", "#4575B4", "#74ADD1",
                                          "#ABD9E9", "#FEE090", "#FDAE61",
                                          "#F46D43", "#D73027", "#A50026"))(100),
                main = paste0("Sample Pairwise Correlation\n",
                             "Min median r = ", round(min(median_cor), 3),
                             " | Flagged: ", sum(cor_outlier)),
                margins = c(8, 8),
                cexRow = max(0.3, min(1, 30 / n_samples)),
                cexCol = max(0.3, min(1, 30 / n_samples)))
    } else {
        # For large sample sets, just show the distribution
        hist(median_cor, breaks = 30,
             main = "Distribution of Median Pairwise Correlations",
             xlab = "Median Spearman correlation",
             col = "steelblue", border = "white")
        abline(v = min_correlation, col = "red", lty = 2, lwd = 2)
        legend("topleft",
               legend = paste0("Threshold = ", min_correlation,
                              "\nFlagged: ", sum(cor_outlier)),
               col = "red", lty = 2, bty = "n")
    }

    dev.off()
    cat("Correlation heatmap written.\n")
}, error = function(e) {
    cat("Correlation heatmap failed:", e$message, "\n")
})

# ---- Generate findings markdown ----
n_any_flag <- sum(results_df$n_flags > 0)
n_multi_flag <- sum(results_df$n_flags >= 2)

findings <- c(
    "# Expression Outlier Detection",
    "",
    "## Parameters",
    "",
    paste0("- PCA SD threshold: ", sd_threshold),
    paste0("- Minimum median correlation: ", min_correlation),
    paste0("- Samples: ", n_samples),
    paste0("- Genes (after filtering): ", nrow(counts_filtered)),
    "",
    "## Summary",
    "",
    paste0("| Detection Method | Outliers | Threshold |"),
    paste0("|-----------------|----------|-----------|"),
    paste0("| PCA distance | ", sum(pca_outlier), " | >", sd_threshold, " SD from centroid |"),
    paste0("| Pairwise correlation | ", sum(cor_outlier), " | median r < ", min_correlation, " |"),
    paste0("| Library size | ", sum(lib_outlier), " | >", sd_threshold, " SD (log10 scale) |"),
    paste0("| **Any flag** | **", n_any_flag, "** | |"),
    paste0("| **Multiple flags** | **", n_multi_flag, "** | |"),
    "",
    "## Library Size Statistics",
    "",
    paste0("- Median: ", sprintf("%.2e", median(lib_sizes))),
    paste0("- Range: ", sprintf("%.2e", min(lib_sizes)), " - ", sprintf("%.2e", max(lib_sizes))),
    paste0("- IQR: ", sprintf("%.2e", quantile(lib_sizes, 0.25)),
           " - ", sprintf("%.2e", quantile(lib_sizes, 0.75))),
    ""
)

# Warnings
warnings_list <- character(0)

if (n_multi_flag > 0) {
    multi_flagged <- results_df[results_df$n_flags >= 2, ]
    warnings_list <- c(warnings_list,
        paste0("- **MULTI-FLAG OUTLIERS**: ", n_multi_flag,
               " sample(s) flagged by multiple detection methods. ",
               "These samples are strong candidates for exclusion or further investigation:"),
        "")
    for (i in seq_len(nrow(multi_flagged))) {
        s <- multi_flagged[i, ]
        flags <- character(0)
        if (s$pca_outlier) flags <- c(flags, "PCA")
        if (s$correlation_outlier) flags <- c(flags, "correlation")
        if (s$library_size_outlier) flags <- c(flags, "library_size")
        warnings_list <- c(warnings_list,
            paste0("  - **", s$sample_id, "**: flagged by ",
                   paste(flags, collapse = ", "),
                   " (median r=", round(s$median_pairwise_correlation, 3),
                   ", PCA z=", round(s$pca_distance_zscore, 2),
                   ", lib_size=", sprintf("%.2e", s$library_size), ")"))
    }
    warnings_list <- c(warnings_list, "")
}

if (n_any_flag > 0 && n_multi_flag == 0) {
    warnings_list <- c(warnings_list,
        paste0("- **SINGLE-FLAG OUTLIERS**: ", n_any_flag,
               " sample(s) flagged by one method. Review these samples but ",
               "single-method flags alone may not warrant exclusion."),
        "")
}

outlier_pct <- n_any_flag / n_samples * 100
if (outlier_pct > 10) {
    warnings_list <- c(warnings_list,
        paste0("- **HIGH OUTLIER RATE**: ", round(outlier_pct, 1),
               "% of samples flagged. This may indicate batch effects, ",
               "systematic technical issues, or overly strict thresholds. ",
               "Consider adjusting thresholds or investigating batch structure."),
        "")
}

if (sum(lib_outlier_low) > 0) {
    low_lib_samples <- results_df$sample_id[results_df$library_size_outlier &
                                              results_df$log10_library_size < lib_mean]
    warnings_list <- c(warnings_list,
        paste0("- **LOW LIBRARY SIZE**: ", length(low_lib_samples),
               " sample(s) have unusually low total counts. ",
               "These may have had poor RNA quality, low input, or sequencing failure. ",
               "Samples: ", paste(low_lib_samples, collapse = ", ")),
        "")
}

if (length(warnings_list) > 0) {
    findings <- c(findings, "## Warnings", "", warnings_list)
} else {
    findings <- c(findings,
        "## Status: PASS",
        "",
        "No expression outliers detected. All samples fall within expected ranges for:",
        "- PCA distance from cohort centroid",
        "- Pairwise correlation with other samples",
        "- Total library size",
        "")
}

# Recommendations
findings <- c(findings,
    "## Recommendations",
    "",
    "- Samples flagged by **multiple methods** are strong candidates for exclusion",
    "- Samples flagged by a single method should be reviewed in context (e.g., known",
    "  biological differences may explain PCA distance)",
    "- If excluding samples, re-run differential expression to assess impact",
    "- Library size outliers may benefit from deeper investigation of RNA quality metrics",
    ""
)

writeLines(findings, file.path(opt$`output-dir`, "expression_outlier_findings.md"))
cat("Findings written to expression_outlier_findings.md\n")

cat("\nExpression outlier detection complete.\n")
