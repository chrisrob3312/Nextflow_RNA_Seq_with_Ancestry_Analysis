#!/usr/bin/env Rscript

# ============================================================================
# Count Normalization - VST, rlog, TMM, TPM
# ============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
    library(optparse)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

option_list <- list(
    make_option("--counts", type = "character", help = "Raw count matrix TSV"),
    make_option("--metadata", type = "character", help = "Sample metadata TSV"),
    make_option("--method", type = "character", default = "vst",
                help = "Normalization method: vst, rlog, tmm, tpm"),
    make_option("--output-prefix", type = "character", default = "normalized_counts"),
    make_option("--plot-dir", type = "character", default = "normalization_qc"),
    make_option("--min-counts", type = "integer", default = 10),
    make_option("--min-samples", type = "integer", default = 3)
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

common <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common]
metadata <- metadata[common, ]

cat("Loaded", nrow(counts), "genes x", ncol(counts), "samples\n")

# ---- Filter low-count genes ----
keep <- rowSums(counts >= opt$`min-counts`) >= opt$`min-samples`
counts_filtered <- counts[keep, ]
cat("After filtering:", nrow(counts_filtered), "genes retained\n")

# ---- Normalize ----
if (opt$method %in% c("vst", "rlog")) {
    dds <- DESeqDataSetFromMatrix(
        countData = round(counts_filtered),
        colData = metadata,
        design = ~ 1
    )
    dds <- estimateSizeFactors(dds)

    if (opt$method == "vst") {
        norm_data <- assay(vst(dds, blind = TRUE))
    } else {
        norm_data <- assay(rlog(dds, blind = TRUE))
    }
} else if (opt$method == "tmm") {
    dge <- DGEList(counts = round(counts_filtered))
    dge <- calcNormFactors(dge, method = "TMM")
    norm_data <- cpm(dge, log = TRUE, prior.count = 1)
}

write.table(norm_data, paste0(opt$`output-prefix`, "_", opt$method, ".tsv"),
            sep = "\t", quote = FALSE)

# ---- Also compute TPM ----
gene_lengths <- rep(1000, nrow(counts_filtered))  # Placeholder; use real lengths if available
rpk <- counts_filtered / gene_lengths * 1000
tpm <- t(t(rpk) / colSums(rpk) * 1e6)
write.table(tpm, paste0(opt$`output-prefix`, "_tpm.tsv"), sep = "\t", quote = FALSE)

# ---- QC plots ----
# PCA
pca <- prcomp(t(norm_data), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:min(5, ncol(pca$x))])
pca_df$sample_id <- rownames(pca_df)
pca_df <- merge(pca_df, metadata, by = "sample_id")
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

for (color_var in c("batch", "cytomolecular_subgroup", "relapse_status", "sex")) {
    if (color_var %in% colnames(pca_df)) {
        p <- ggplot(pca_df, aes_string(x = "PC1", y = "PC2", color = color_var)) +
            geom_point(size = 3, alpha = 0.8) +
            labs(x = paste0("PC1 (", var_explained[1], "%)"),
                 y = paste0("PC2 (", var_explained[2], "%)"),
                 title = paste("PCA -", color_var)) +
            theme_bw() +
            theme(legend.position = "right")
        ggsave(file.path(opt$`plot-dir`, paste0("pca_", color_var, ".pdf")),
               p, width = 8, height = 6)
    }
}

# Sample correlation heatmap
cor_mat <- cor(norm_data, method = "spearman")
annot_df <- metadata[, intersect(c("batch", "cytomolecular_subgroup", "relapse_status"), colnames(metadata)), drop = FALSE]
pdf(file.path(opt$`plot-dir`, "sample_correlation_heatmap.pdf"), width = 12, height = 10)
pheatmap(cor_mat,
         annotation_col = annot_df,
         color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         main = "Sample Correlation (Spearman)")
dev.off()

# Density plot
pdf(file.path(opt$`plot-dir`, "expression_density.pdf"), width = 10, height = 6)
plot(density(norm_data[, 1]), main = "Expression Density", xlab = "Normalized Expression",
     ylim = c(0, 0.3), col = rgb(0, 0, 1, 0.3))
for (i in 2:ncol(norm_data)) {
    lines(density(norm_data[, i]), col = rgb(0, 0, 1, 0.3))
}
dev.off()

cat("Normalization complete. Output written with prefix:", opt$`output-prefix`, "\n")
