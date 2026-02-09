#!/usr/bin/env Rscript

# ============================================================================
# Molecular Subtyping / Cell-of-Origin Classification
# For B-ALL: Uses gene expression signatures to classify molecular subtypes
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(pheatmap)
    library(class)
    library(randomForest)
})

option_list <- list(
    make_option("--expression", type = "character"),
    make_option("--metadata", type = "character"),
    make_option("--method", type = "character", default = "consensus"),
    make_option("--output-dir", type = "character", default = "subtyping_results"),
    make_option("--plot-dir", type = "character", default = "subtyping_plots")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)

expr <- read.delim(opt$expression, row.names = 1, check.names = FALSE)
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

common <- intersect(colnames(expr), rownames(metadata))
expr <- expr[, common]
metadata <- metadata[common, ]

# ---- B-ALL specific molecular subtyping gene signatures ----
# Key B-ALL subtypes and their marker genes
ball_signatures <- list(
    "ETV6-RUNX1" = c("ETV6", "RUNX1", "CLIC5", "AGAP1", "AUTS2"),
    "BCR-ABL1" = c("BCR", "ABL1", "SEMA6A", "BMP6", "PRG2"),
    "KMT2A" = c("KMT2A", "MEIS1", "HOXA9", "HOXA10", "FLT3"),
    "Hyperdiploid" = c("CRLF2", "TSLP", "IL3RA"),
    "TCF3-PBX1" = c("TCF3", "PBX1", "WNT16", "CTNNAL1"),
    "iAMP21" = c("RUNX1", "ERG", "DYRK1A"),
    "Ph-like" = c("CRLF2", "JAK2", "EPOR", "ABL1", "ABL2", "CSF1R", "PDGFRB"),
    "DUX4" = c("DUX4", "ERG", "AGAP1"),
    "PAX5alt" = c("PAX5", "DACH1", "ID4"),
    "MEF2D" = c("MEF2D", "BCL9", "HDAC9"),
    "ZNF384" = c("ZNF384", "CLCF1", "BMP6"),
    "NUTM1" = c("NUTM1", "BRD4")
)

# Score each sample for each subtype signature
scores <- data.frame(sample_id = common)
for (subtype in names(ball_signatures)) {
    genes <- ball_signatures[[subtype]]
    found_genes <- genes[genes %in% rownames(expr)]
    if (length(found_genes) > 0) {
        scores[[subtype]] <- colMeans(expr[found_genes, , drop = FALSE], na.rm = TRUE)
    } else {
        scores[[subtype]] <- NA
    }
}

# Classify based on highest signature score
score_cols <- setdiff(colnames(scores), "sample_id")
if (length(score_cols) > 0) {
    score_matrix <- as.matrix(scores[, score_cols])
    scores$predicted_subtype <- score_cols[apply(score_matrix, 1, which.max)]
    scores$max_score <- apply(score_matrix, 1, max, na.rm = TRUE)
    scores$confidence <- apply(score_matrix, 1, function(x) {
        sorted <- sort(x, decreasing = TRUE, na.last = TRUE)
        if (length(sorted) >= 2 && !is.na(sorted[1]) && !is.na(sorted[2]) && sorted[2] != 0) {
            sorted[1] / sorted[2]
        } else {
            NA
        }
    })
}

write.table(scores, file.path(opt$`output-dir`, "subtypes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(scores[, c("sample_id", score_cols)],
            file.path(opt$`output-dir`, "classifier_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Plots ----
# Heatmap of subtype scores
if (length(score_cols) > 0) {
    plot_mat <- t(as.matrix(scores[, score_cols]))
    colnames(plot_mat) <- scores$sample_id

    annot_df <- metadata[scores$sample_id, intersect(c("cytomolecular_subgroup", "relapse_status"), colnames(metadata)), drop = FALSE]

    pdf(file.path(opt$`plot-dir`, "subtype_scores_heatmap.pdf"), width = 14, height = 8)
    pheatmap(plot_mat, scale = "row",
             annotation_col = annot_df,
             main = "B-ALL Molecular Subtype Signature Scores",
             fontsize_col = 7, fontsize_row = 9)
    dev.off()

    # Compare predicted vs known
    if ("cytomolecular_subgroup" %in% colnames(metadata)) {
        comparison <- data.frame(
            sample_id = scores$sample_id,
            known = metadata[scores$sample_id, "cytomolecular_subgroup"],
            predicted = scores$predicted_subtype,
            confidence = scores$confidence
        )
        write.table(comparison, file.path(opt$`output-dir`, "subtype_comparison.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

cat("Molecular subtyping complete.\n")
