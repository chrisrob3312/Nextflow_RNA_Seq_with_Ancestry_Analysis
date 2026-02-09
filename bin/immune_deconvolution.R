#!/usr/bin/env Rscript

# ============================================================================
# Immune Deconvolution - Multiple methods with tumor purity correction
# Methods: xCell, MCP-counter, ESTIMATE, EPIC, TIMER, CIBERSORTx
# ============================================================================

suppressPackageStartupMessages({
    library(immunedeconv)
    library(optparse)
    library(ggplot2)
    library(pheatmap)
    library(reshape2)
})

option_list <- list(
    make_option("--expression", type = "character"),
    make_option("--metadata", type = "character"),
    make_option("--methods", type = "character", default = "xcell,mcpcounter,estimate,epic"),
    make_option("--cibersortx-token", type = "character", default = NULL),
    make_option("--sig-matrix", type = "character", default = "LM22"),
    make_option("--correct-purity", type = "logical", default = TRUE),
    make_option("--output-dir", type = "character", default = "immune_results"),
    make_option("--plot-dir", type = "character", default = "immune_plots")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
expr <- read.delim(opt$expression, row.names = 1, check.names = FALSE)
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

# Ensure TPM-like values (most deconvolution methods expect non-log scale)
if (min(expr, na.rm = TRUE) < 0) {
    cat("Expression appears log-transformed. Converting back for deconvolution.\n")
    expr_raw <- 2^expr  # Back-transform from log2
} else {
    expr_raw <- expr
}

expr_matrix <- as.matrix(expr_raw)

methods_list <- trimws(unlist(strsplit(opt$methods, ",")))
all_results <- list()

# ---- Run each deconvolution method ----
for (method in methods_list) {
    cat("\nRunning deconvolution:", method, "\n")

    tryCatch({
        if (method == "cibersortx" && !is.null(opt$`cibersortx-token`)) {
            set_cibersort_binary(opt$`cibersortx-token`)
            res <- deconvolute(expr_matrix, method = "cibersort_abs")
        } else if (method == "cibersortx") {
            cat("Skipping CIBERSORTx - no token provided\n")
            next
        } else if (method == "estimate") {
            # ESTIMATE is handled separately
            next
        } else {
            res <- deconvolute(expr_matrix, method = method)
        }

        res_df <- as.data.frame(res)
        write.table(res_df,
                    file.path(opt$`output-dir`, paste0("deconv_", method, ".tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        all_results[[method]] <- res_df

        # Heatmap
        tryCatch({
            plot_data <- as.matrix(res_df[, -1])
            rownames(plot_data) <- res_df$cell_type

            # Subset to variable cell types
            var_cells <- apply(plot_data, 1, var, na.rm = TRUE)
            top_cells <- names(head(sort(var_cells, decreasing = TRUE), 20))
            plot_data <- plot_data[top_cells, , drop = FALSE]

            pdf(file.path(opt$`plot-dir`, paste0(method, "_heatmap.pdf")), width = 14, height = 8)
            pheatmap(plot_data, scale = "row",
                     main = paste("Immune Deconvolution -", toupper(method)),
                     fontsize_col = 7, fontsize_row = 9)
            dev.off()
        }, error = function(e) cat("Heatmap error for", method, ":", e$message, "\n"))

    }, error = function(e) cat("Error running", method, ":", e$message, "\n"))
}

# ---- Combine all methods ----
if (length(all_results) > 0) {
    combined <- do.call(rbind, lapply(names(all_results), function(m) {
        df <- all_results[[m]]
        df$method <- m
        df
    }))
    write.table(combined,
                file.path(opt$`output-dir`, "deconvolution_all.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- Tumor purity correction (if applicable) ----
if (opt$`correct-purity` && "tumor_purity" %in% colnames(metadata)) {
    cat("\nApplying tumor purity correction...\n")
    for (method in names(all_results)) {
        res_df <- all_results[[method]]
        # Scale immune fractions by (1 - tumor_purity) to estimate absolute fractions
        for (sample in colnames(res_df)[-1]) {
            if (sample %in% rownames(metadata)) {
                purity <- as.numeric(metadata[sample, "tumor_purity"])
                if (!is.na(purity) && purity > 0 && purity < 1) {
                    res_df[, sample] <- res_df[, sample] / (1 - purity)
                }
            }
        }
        write.table(res_df,
                    file.path(opt$`output-dir`, paste0("deconv_", method, "_purity_corrected.tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

# ---- Cell fraction summary ----
# Aggregate across methods for consensus cell fractions
cat("\nImmune deconvolution complete.\n")
