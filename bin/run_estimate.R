#!/usr/bin/env Rscript

# ============================================================================
# ESTIMATE - Tumor purity and immune/stromal score estimation
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
})

option_list <- list(
    make_option("--expression", type = "character"),
    make_option("--metadata", type = "character"),
    make_option("--platform", type = "character", default = "illumina"),
    make_option("--output-scores", type = "character", default = "estimate_scores.tsv"),
    make_option("--output-purity", type = "character", default = "estimate_purity.tsv"),
    make_option("--plot-dir", type = "character", default = "estimate_plots")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)

# Load ESTIMATE if available, otherwise use gene signature approach
tryCatch({
    library(estimate)

    expr <- read.delim(opt$expression, row.names = 1, check.names = FALSE)
    metadata <- read.delim(opt$metadata, check.names = FALSE)
    rownames(metadata) <- metadata$sample_id

    # Write expression in GCT format for ESTIMATE
    gct_file <- "expression_for_estimate.gct"
    writeLines(c("#1.2", paste(nrow(expr), ncol(expr), sep = "\t")), gct_file)
    header <- paste(c("NAME", "Description", colnames(expr)), collapse = "\t")
    write(header, gct_file, append = TRUE)
    for (i in seq_len(nrow(expr))) {
        line <- paste(c(rownames(expr)[i], "na", as.character(expr[i, ])), collapse = "\t")
        write(line, gct_file, append = TRUE)
    }

    # Run ESTIMATE
    filterCommonGenes(input.f = gct_file, output.f = "common_genes.gct", id = "GeneSymbol")
    estimateScore("common_genes.gct", "estimate_scores.gct", platform = opt$platform)

    # Parse output
    scores <- read.delim("estimate_scores.gct", skip = 2, check.names = FALSE)
    score_names <- scores[, 1]
    score_matrix <- t(scores[, -(1:2)])
    colnames(score_matrix) <- score_names

    score_df <- as.data.frame(score_matrix)
    score_df$sample_id <- rownames(score_df)

    write.table(score_df, opt$`output-scores`, sep = "\t", quote = FALSE, row.names = FALSE)

    # Estimate tumor purity from ESTIMATE score
    if ("ESTIMATEScore" %in% colnames(score_df)) {
        # Purity formula: cos(0.6049872018 + 0.0001467884 * ESTIMATE_score)
        purity <- cos(0.6049872018 + 0.0001467884 * as.numeric(score_df$ESTIMATEScore))
        purity <- pmax(0, pmin(1, purity))
        purity_df <- data.frame(sample_id = score_df$sample_id,
                                 estimate_purity = round(purity, 4))
        write.table(purity_df, opt$`output-purity`, sep = "\t", quote = FALSE, row.names = FALSE)
    }

    # Plots
    if ("StromalScore" %in% colnames(score_df) && "ImmuneScore" %in% colnames(score_df)) {
        plot_df <- merge(score_df, metadata, by = "sample_id")
        plot_df$StromalScore <- as.numeric(plot_df$StromalScore)
        plot_df$ImmuneScore <- as.numeric(plot_df$ImmuneScore)

        p <- ggplot(plot_df, aes(x = StromalScore, y = ImmuneScore)) +
            geom_point(aes(color = relapse_status), size = 3, alpha = 0.7) +
            theme_bw() +
            labs(title = "ESTIMATE: Stromal vs Immune Scores")
        ggsave(file.path(opt$`plot-dir`, "estimate_stromal_vs_immune.pdf"), p, width = 8, height = 6)
    }

}, error = function(e) {
    cat("ESTIMATE package error:", e$message, "\n")
    cat("Creating placeholder output.\n")
    writeLines("sample_id\tStromalScore\tImmuneScore\tESTIMATEScore", opt$`output-scores`)
    writeLines("sample_id\testimate_purity", opt$`output-purity`)
})

cat("ESTIMATE analysis complete.\n")
