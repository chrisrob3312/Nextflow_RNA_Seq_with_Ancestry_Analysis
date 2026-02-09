#!/usr/bin/env Rscript

# ============================================================================
# Batch Effect Correction - ComBat-seq (for counts) and SVA
# ============================================================================

suppressPackageStartupMessages({
    library(sva)
    library(optparse)
    library(ggplot2)
    library(DESeq2)
})

option_list <- list(
    make_option("--counts", type = "character"),
    make_option("--metadata", type = "character"),
    make_option("--method", type = "character", default = "combat_seq"),
    make_option("--batch-variable", type = "character", default = "batch"),
    make_option("--output", type = "character", default = "batch_corrected_counts.tsv"),
    make_option("--plot-dir", type = "character", default = "batch_correction_qc")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

common <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common]
metadata <- metadata[common, ]

batch <- metadata[[opt$`batch-variable`]]

if (length(unique(batch[!is.na(batch)])) <= 1) {
    cat("Only one batch detected. Skipping batch correction.\n")
    write.table(counts, opt$output, sep = "\t", quote = FALSE)
    quit(save = "no")
}

cat("Batch variable:", opt$`batch-variable`, "\n")
cat("Batches:", paste(unique(batch), collapse = ", "), "\n")

# PCA before correction
dds_pre <- DESeqDataSetFromMatrix(round(counts), metadata, design = ~ 1)
vsd_pre <- assay(vst(dds_pre, blind = TRUE))
pca_pre <- prcomp(t(vsd_pre), scale. = TRUE)

if (opt$method == "combat_seq") {
    cat("Running ComBat-seq...\n")
    corrected <- ComBat_seq(as.matrix(round(counts)), batch = batch,
                             group = NULL, covar_mod = NULL)
} else if (opt$method == "sva") {
    cat("Running SVA to estimate surrogate variables...\n")
    mod <- model.matrix(~ 1, data = metadata)
    svobj <- sva(as.matrix(vsd_pre), mod, method = "irw")
    cat("Estimated", svobj$n.sv, "surrogate variables\n")
    # For SVA, we output the surrogate variables to be used as covariates
    sv_df <- data.frame(sample_id = rownames(metadata))
    for (i in seq_len(svobj$n.sv)) {
        sv_df[[paste0("SV", i)]] <- svobj$sv[, i]
    }
    write.table(sv_df, gsub("\\.tsv$", "_surrogate_variables.tsv", opt$output),
                sep = "\t", quote = FALSE, row.names = FALSE)
    corrected <- counts  # SVA doesn't directly correct counts
}

write.table(corrected, opt$output, sep = "\t", quote = FALSE)

# PCA after correction
dds_post <- DESeqDataSetFromMatrix(round(as.matrix(corrected)), metadata, design = ~ 1)
vsd_post <- assay(vst(dds_post, blind = TRUE))
pca_post <- prcomp(t(vsd_post), scale. = TRUE)

# Plot comparison
pdf(file.path(opt$`plot-dir`, "batch_correction_pca.pdf"), width = 14, height = 6)
par(mfrow = c(1, 2))
cols <- as.numeric(as.factor(batch))
plot(pca_pre$x[, 1], pca_pre$x[, 2], col = cols, pch = 19,
     main = "Before Correction", xlab = "PC1", ylab = "PC2")
legend("topright", legend = unique(batch), col = unique(cols), pch = 19, cex = 0.8)
plot(pca_post$x[, 1], pca_post$x[, 2], col = cols, pch = 19,
     main = "After Correction", xlab = "PC1", ylab = "PC2")
legend("topright", legend = unique(batch), col = unique(cols), pch = 19, cex = 0.8)
dev.off()

cat("Batch correction complete.\n")
