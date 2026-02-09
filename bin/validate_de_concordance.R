#!/usr/bin/env Rscript

# ============================================================================
# DESeq2 vs limma-voom Concordance Validation
# Compares differential expression results from both methods to assess
# reproducibility and flag potential method-specific artifacts.
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
})

# ---- Parse arguments ----
option_list <- list(
    make_option("--deseq2", type = "character",
                help = "DESeq2 results TSV file (gene rownames, must have log2FoldChange and padj columns)"),
    make_option("--limma", type = "character",
                help = "limma-voom results TSV file (gene rownames, must have logFC and adj.P.Val columns)"),
    make_option("--contrast", type = "character", default = "unknown",
                help = "Contrast name for labeling outputs [default: %default]"),
    make_option("--output-dir", type = "character", default = "concordance_results",
                help = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate inputs ----
if (is.null(opt$deseq2) || is.null(opt$limma)) {
    stop("Both --deseq2 and --limma arguments are required.")
}

if (!file.exists(opt$deseq2)) {
    stop(paste("DESeq2 results file not found:", opt$deseq2))
}

if (!file.exists(opt$limma)) {
    stop(paste("limma-voom results file not found:", opt$limma))
}

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- Load results ----
cat("Loading DESeq2 results:", opt$deseq2, "\n")
deseq2_res <- read.delim(opt$deseq2, check.names = FALSE, row.names = 1)

cat("Loading limma-voom results:", opt$limma, "\n")
limma_res <- read.delim(opt$limma, check.names = FALSE, row.names = 1)

# Validate that required columns exist
deseq2_lfc_col <- NULL
deseq2_padj_col <- NULL
limma_lfc_col <- NULL
limma_padj_col <- NULL

if ("log2FoldChange" %in% colnames(deseq2_res)) {
    deseq2_lfc_col <- "log2FoldChange"
} else {
    stop("DESeq2 results missing 'log2FoldChange' column.")
}

if ("padj" %in% colnames(deseq2_res)) {
    deseq2_padj_col <- "padj"
} else if ("adj.P.Val" %in% colnames(deseq2_res)) {
    deseq2_padj_col <- "adj.P.Val"
} else {
    stop("DESeq2 results missing 'padj' or 'adj.P.Val' column.")
}

if ("logFC" %in% colnames(limma_res)) {
    limma_lfc_col <- "logFC"
} else if ("log2FoldChange" %in% colnames(limma_res)) {
    limma_lfc_col <- "log2FoldChange"
} else {
    stop("limma results missing 'logFC' or 'log2FoldChange' column.")
}

if ("adj.P.Val" %in% colnames(limma_res)) {
    limma_padj_col <- "adj.P.Val"
} else if ("padj" %in% colnames(limma_res)) {
    limma_padj_col <- "padj"
} else {
    stop("limma results missing 'adj.P.Val' or 'padj' column.")
}

# ---- Merge by gene ----
common_genes <- intersect(rownames(deseq2_res), rownames(limma_res))
cat("Common genes between methods:", length(common_genes), "\n")
cat("  DESeq2 genes:", nrow(deseq2_res), "\n")
cat("  limma genes:", nrow(limma_res), "\n")

if (length(common_genes) < 10) {
    warning("Fewer than 10 common genes found. Concordance metrics may be unreliable.")
}

if (length(common_genes) == 0) {
    # Write empty results and warning findings
    findings <- c(
        paste0("# DE Concordance Validation: ", opt$contrast),
        "",
        "## Status: FAILED",
        "",
        "**No common genes found between DESeq2 and limma-voom results.**",
        "",
        "This may indicate:",
        "- Different gene ID formats between the two result files",
        "- One or both analyses produced empty results",
        "- Input files may be corrupted or malformed",
        ""
    )
    writeLines(findings, file.path(opt$`output-dir`, "concordance_findings.md"))
    cat("No common genes. Findings written.\n")
    quit(status = 0)
}

merged <- data.frame(
    gene = common_genes,
    deseq2_lfc = deseq2_res[common_genes, deseq2_lfc_col],
    limma_lfc = limma_res[common_genes, limma_lfc_col],
    deseq2_padj = deseq2_res[common_genes, deseq2_padj_col],
    limma_padj = limma_res[common_genes, limma_padj_col],
    stringsAsFactors = FALSE
)

# Remove rows with NA in fold changes
merged_complete <- merged[!is.na(merged$deseq2_lfc) & !is.na(merged$limma_lfc), ]
cat("Genes with complete fold-change data:", nrow(merged_complete), "\n")

# ---- Compute concordance metrics ----

# 1. Spearman correlation of log2FC
spearman_cor <- cor(merged_complete$deseq2_lfc, merged_complete$limma_lfc,
                    method = "spearman", use = "complete.obs")
pearson_cor <- cor(merged_complete$deseq2_lfc, merged_complete$limma_lfc,
                   method = "pearson", use = "complete.obs")

cat("Spearman correlation of log2FC:", round(spearman_cor, 4), "\n")
cat("Pearson correlation of log2FC:", round(pearson_cor, 4), "\n")

# 2. Jaccard index of significant genes (padj < 0.05)
sig_deseq2 <- merged$gene[!is.na(merged$deseq2_padj) & merged$deseq2_padj < 0.05]
sig_limma <- merged$gene[!is.na(merged$limma_padj) & merged$limma_padj < 0.05]

n_intersection <- length(intersect(sig_deseq2, sig_limma))
n_union <- length(union(sig_deseq2, sig_limma))

jaccard <- if (n_union > 0) n_intersection / n_union else NA
cat("Significant genes - DESeq2:", length(sig_deseq2),
    "| limma:", length(sig_limma),
    "| overlap:", n_intersection, "\n")
cat("Jaccard index:", round(jaccard, 4), "\n")

# 3. Direction agreement among jointly significant genes
joint_sig <- intersect(sig_deseq2, sig_limma)
if (length(joint_sig) > 0) {
    joint_data <- merged[merged$gene %in% joint_sig, ]
    same_direction <- sum(sign(joint_data$deseq2_lfc) == sign(joint_data$limma_lfc), na.rm = TRUE)
    direction_agreement <- same_direction / length(joint_sig) * 100
} else {
    same_direction <- 0
    direction_agreement <- NA
}
cat("Direction agreement (jointly significant):", round(direction_agreement, 2), "%\n")

# 4. Direction agreement among all genes with non-zero fold changes
nonzero <- merged_complete[merged_complete$deseq2_lfc != 0 & merged_complete$limma_lfc != 0, ]
if (nrow(nonzero) > 0) {
    all_direction_agreement <- sum(sign(nonzero$deseq2_lfc) == sign(nonzero$limma_lfc)) /
                                nrow(nonzero) * 100
} else {
    all_direction_agreement <- NA
}

# ---- Write concordance summary ----
summary_df <- data.frame(
    contrast = opt$contrast,
    n_common_genes = length(common_genes),
    n_deseq2_only = nrow(deseq2_res) - length(common_genes),
    n_limma_only = nrow(limma_res) - length(common_genes),
    spearman_correlation = round(spearman_cor, 4),
    pearson_correlation = round(pearson_cor, 4),
    n_sig_deseq2 = length(sig_deseq2),
    n_sig_limma = length(sig_limma),
    n_sig_both = n_intersection,
    jaccard_index = round(jaccard, 4),
    n_joint_sig = length(joint_sig),
    n_same_direction = same_direction,
    direction_agreement_pct = round(direction_agreement, 2),
    all_genes_direction_agreement_pct = round(all_direction_agreement, 2),
    stringsAsFactors = FALSE
)

write.table(summary_df, file.path(opt$`output-dir`, "concordance_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary written to concordance_summary.tsv\n")

# ---- Generate scatterplot ----
tryCatch({
    pdf(file.path(opt$`output-dir`, "concordance_scatterplot.pdf"), width = 8, height = 8)

    # Classify genes for coloring
    merged_complete$category <- "Not significant"
    merged_complete$category[merged_complete$gene %in% sig_deseq2 &
                              !(merged_complete$gene %in% sig_limma)] <- "DESeq2 only"
    merged_complete$category[!(merged_complete$gene %in% sig_deseq2) &
                              merged_complete$gene %in% sig_limma] <- "limma only"
    merged_complete$category[merged_complete$gene %in% joint_sig] <- "Both significant"
    merged_complete$category <- factor(merged_complete$category,
                                        levels = c("Not significant", "DESeq2 only",
                                                   "limma only", "Both significant"))

    p <- ggplot(merged_complete, aes(x = deseq2_lfc, y = limma_lfc, color = category)) +
        geom_point(alpha = 0.4, size = 0.8) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
        scale_color_manual(values = c("Not significant" = "grey70",
                                       "DESeq2 only" = "#E69F00",
                                       "limma only" = "#56B4E9",
                                       "Both significant" = "#D55E00")) +
        labs(
            title = paste0("DESeq2 vs limma-voom: ", opt$contrast),
            subtitle = paste0("Spearman rho = ", round(spearman_cor, 3),
                             " | Jaccard = ", round(jaccard, 3),
                             " | Direction agreement = ", round(direction_agreement, 1), "%"),
            x = "DESeq2 log2FoldChange",
            y = "limma-voom logFC",
            color = "Significance"
        ) +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom")

    print(p)
    dev.off()
    cat("Scatterplot written to concordance_scatterplot.pdf\n")
}, error = function(e) {
    cat("Scatterplot generation failed:", e$message, "\n")
})

# ---- Generate findings markdown ----
findings <- c(
    paste0("# DE Concordance Validation: ", opt$contrast),
    "",
    "## Summary Metrics",
    "",
    paste0("| Metric | Value |"),
    paste0("|--------|-------|"),
    paste0("| Common genes | ", length(common_genes), " |"),
    paste0("| Spearman correlation (log2FC) | ", round(spearman_cor, 4), " |"),
    paste0("| Pearson correlation (log2FC) | ", round(pearson_cor, 4), " |"),
    paste0("| Significant genes (DESeq2, padj<0.05) | ", length(sig_deseq2), " |"),
    paste0("| Significant genes (limma, padj<0.05) | ", length(sig_limma), " |"),
    paste0("| Significant in both | ", n_intersection, " |"),
    paste0("| Jaccard index | ", round(jaccard, 4), " |"),
    paste0("| Direction agreement (jointly sig.) | ",
           ifelse(is.na(direction_agreement), "N/A", paste0(round(direction_agreement, 1), "%")), " |"),
    paste0("| Direction agreement (all genes) | ",
           ifelse(is.na(all_direction_agreement), "N/A", paste0(round(all_direction_agreement, 1), "%")), " |"),
    ""
)

# Warnings section
warnings_list <- character(0)

if (!is.na(spearman_cor) && spearman_cor < 0.7) {
    warnings_list <- c(warnings_list,
        paste0("- **LOW CORRELATION**: Spearman correlation (", round(spearman_cor, 3),
               ") is below 0.7. The two methods show substantial disagreement in ",
               "fold-change estimates. This may indicate sensitivity to normalization ",
               "or model assumptions. Investigate genes with large discrepancies."))
}

if (!is.na(jaccard) && jaccard < 0.5) {
    warnings_list <- c(warnings_list,
        paste0("- **LOW OVERLAP**: Jaccard index (", round(jaccard, 3),
               ") is below 0.5. Fewer than half of significant genes are shared. ",
               "Consider focusing downstream analyses on the concordant gene set ",
               "(n=", n_intersection, ") for more robust conclusions."))
}

if (!is.na(direction_agreement) && direction_agreement < 95) {
    warnings_list <- c(warnings_list,
        paste0("- **DIRECTION DISAGREEMENT**: ", round(100 - direction_agreement, 1),
               "% of jointly significant genes have opposite fold-change direction. ",
               "These discordant genes should be examined individually."))
}

if (length(common_genes) < 100) {
    warnings_list <- c(warnings_list,
        paste0("- **FEW COMMON GENES**: Only ", length(common_genes),
               " genes were tested by both methods. Results may be unreliable. ",
               "Check that gene ID formats match between inputs."))
}

if (length(sig_deseq2) == 0 && length(sig_limma) == 0) {
    warnings_list <- c(warnings_list,
        "- **NO SIGNIFICANT GENES**: Neither method found significant genes at padj<0.05. ",
        "  This may indicate insufficient power or a truly null contrast.")
}

if (length(warnings_list) > 0) {
    findings <- c(findings, "## Warnings", "", warnings_list, "")
} else {
    findings <- c(findings,
        "## Status: PASS",
        "",
        "Concordance between DESeq2 and limma-voom is within acceptable thresholds.",
        "Both methods show good agreement in fold-change estimates and significant gene sets.",
        "")
}

# Interpretation
findings <- c(findings,
    "## Interpretation",
    "",
    "- **Spearman correlation >= 0.7**: Methods agree well on relative ranking of fold changes.",
    "- **Jaccard index >= 0.5**: Majority of significant genes are detected by both methods.",
    "- **Direction agreement >= 95%**: Jointly significant genes show consistent up/down regulation.",
    "",
    "Concordant genes (significant in both methods with same direction) represent the most",
    "robust differential expression signals and should be prioritized for biological interpretation.",
    ""
)

writeLines(findings, file.path(opt$`output-dir`, "concordance_findings.md"))
cat("Findings written to concordance_findings.md\n")

cat("\nDE concordance validation complete.\n")
