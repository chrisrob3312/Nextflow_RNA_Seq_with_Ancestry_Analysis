#!/usr/bin/env Rscript
# ============================================================================
# 18 - Molecular Subtyping (B-ALL)
# ============================================================================
# Classifies B-ALL samples using gene signature scoring (GSVA/ssGSEA).
# Defines subtype signatures for: ETV6-RUNX1, BCR-ABL1, KMT2A, Hyperdiploid,
# Ph-like, DUX4, iAMP21, PAX5alt, MEF2D, ZNF384, NUTM1, and others.
# Assigns subtype by highest scoring signature above confidence threshold.
#
# Required packages: GSVA, GSEABase, ggplot2, pheatmap, tidyverse, optparse
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(GSVA)
    library(tidyverse)
    library(ggplot2)
    library(pheatmap)
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
    make_option("--confidence_threshold", type = "numeric", default = 0.3,
                help = "Minimum score to assign a subtype [default: %default]"),
    make_option("--threads", type = "integer", default = 4,
                help = "Number of threads [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$normalized_counts) || is.null(opt$metadata) || is.null(opt$output_dir)) {
    stop("Required arguments: --normalized_counts, --metadata, --output_dir")
}

cat("=== B-ALL Molecular Subtyping ===\n")
cat("Normalized counts:", opt$normalized_counts, "\n")
cat("Metadata:", opt$metadata, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Confidence threshold:", opt$confidence_threshold, "\n")
cat("Threads:", opt$threads, "\n\n")

# ============================================================================
# Define B-ALL Subtype Gene Signatures
# ============================================================================
cat("Defining B-ALL subtype gene signatures...\n")

ball_signatures <- list(
    # ETV6-RUNX1 (t(12;21)) - most common in pediatric B-ALL
    ETV6_RUNX1 = c("ETV6", "RUNX1", "CLIC5", "TPMT", "PDE4B", "NR3C1",
                    "ARHGAP24", "EPOR", "IGFBP7", "ATP8B4", "PTPRK",
                    "CACNA2D1", "ASNS", "C1orf116", "RAG1", "DNTT",
                    "CD44", "PBX1", "TCL1A", "DDIT4L", "MDFIC",
                    "RPS6KA2", "AIM2", "KCNK12"),

    # BCR-ABL1 (Ph+) - t(9;22)
    BCR_ABL1 = c("BCR", "ABL1", "SEMA6A", "PON2", "IGSF1", "BLK",
                  "IFITM1", "TLE4", "MYB", "TSPAN7", "IL2RA",
                  "CDKN2A", "CD9", "SCARB1", "SLC2A5", "ORAI2",
                  "MEIS1", "BMPR1B", "GATA3", "CHN2", "DUSP6",
                  "SPRY2", "KCNMA1", "S100A6"),

    # KMT2A (MLL) rearranged - t(4;11), t(9;11), t(11;19)
    KMT2A = c("KMT2A", "HOXA9", "HOXA10", "HOXA7", "HOXA5", "MEIS1",
              "PBX3", "FLT3", "LAMP5", "CPNE8", "IRX1", "ZEB1",
              "LGALS1", "CCNA1", "ADAM10", "ENPP1", "PROSER2",
              "NEGR1", "CDK6", "SENP6", "BEX1", "MYO1E",
              "PROM1", "LMO4"),

    # Hyperdiploid (51-67 chromosomes)
    Hyperdiploid = c("STMN1", "GINS2", "TOP2A", "CENPF", "CRLF2",
                     "GPR56", "CD44", "HMGA2", "MDFIC", "FHL1",
                     "MUC4", "IL3RA", "EPOR", "TSPAN7", "KCNJ15",
                     "ID4", "NOTCH1", "PTPRM", "GPR171", "BMPR1B",
                     "ADARB1", "SDK1", "CLCA4", "SYT1"),

    # Ph-like (BCR-ABL1-like, CRLF2/JAK pathway)
    Ph_like = c("CRLF2", "JAK2", "TSLP", "IL7R", "SH2B3", "EPOR",
                "ABL1", "ABL2", "CSF1R", "PDGFRB", "NTRK3", "FLT3",
                "IL2RB", "BLNK", "VPREB1", "ZCCHC7", "SOCS2",
                "CA6", "IFITM1", "MUC4", "IGJ", "PON2",
                "CHN2", "BMPR1B"),

    # DUX4 rearranged (IGH-DUX4)
    DUX4 = c("DUX4", "DUXAP8", "AGAP1", "ERG", "ETS2", "NFATC4",
             "CD2", "CLEC12A", "RASGRP1", "SOCS2", "GATA3",
             "CLIC5", "VIM", "LPXN", "DDIT4L", "HLA-DQB1",
             "ALDOC", "ZEB2", "PYHIN1", "NRXN3", "CDH2",
             "RGS1", "FHIT", "TMEM156"),

    # iAMP21 (intrachromosomal amplification of chromosome 21)
    iAMP21 = c("RUNX1", "ERG", "ETS2", "DYRK1A", "CHAF1B", "SON",
               "HMGN1", "BRWD1", "BACE2", "SIM2", "HLCS",
               "FTCD", "CSTB", "PRMT2", "RRP1", "OLIG1",
               "OLIG2", "BACH1", "TIAM1", "GABPA", "PKNOX1",
               "CBS", "SLC19A1", "SOD1"),

    # PAX5 alterations (PAX5alt - fusions/mutations)
    PAX5alt = c("PAX5", "EBF1", "BACH2", "ID4", "SOX4", "FOXO1",
                "VPREB1", "CD79A", "CD79B", "BLK", "BLNK",
                "IKZF1", "IKZF3", "IRF4", "PRDM1", "TCF3",
                "LEF1", "ROR1", "CD19", "MS4A1", "BANK1",
                "RAG1", "RAG2", "IGHM"),

    # MEF2D rearranged
    MEF2D = c("MEF2D", "HDAC9", "GRIA4", "NRXN3", "STMN1",
              "BCL9", "SS18", "DAZAP1", "HNRNPUL1", "FOXP1",
              "PBX1", "ITM2C", "NEGR1", "SLC1A4", "CLIC5",
              "GATA3", "MYB", "VPREB1", "IRF4", "MEIS2",
              "CDK6", "ROR1", "IGLL1", "FBXO16"),

    # ZNF384 rearranged (EP300-ZNF384, TCF3-ZNF384, etc.)
    ZNF384 = c("ZNF384", "BMP6", "CLCF1", "VIM", "NRXN3", "CHI3L1",
               "RUNX2", "CEBPA", "CEBPB", "CD13", "CD33",
               "ITGAM", "LYZ", "MPO", "CSF3R", "ANPEP",
               "FLT3", "WT1", "GATA2", "GRIA4", "STMN1",
               "FOXP1", "SOX4", "PROM1"),

    # NUTM1 rearranged
    NUTM1 = c("NUTM1", "BRD4", "BRD3", "NSD3", "ACLY", "SLC1A3",
              "CDK6", "HOXA9", "HOXA7", "MEIS1", "PBX3",
              "HMGA2", "SOX11", "CD34", "FLT3", "KIT",
              "GATA2", "ERG", "SPI1", "CBFA2T3", "NPM1",
              "MYC", "BCL2", "MCL1")
)

cat(sprintf("  Defined %d subtype signatures\n", length(ball_signatures)))
for (sig in names(ball_signatures)) {
    cat(sprintf("    %s: %d genes\n", sig, length(ball_signatures[[sig]])))
}

# ============================================================================
# Load data
# ============================================================================
cat("\nLoading normalized count matrix...\n")
counts <- read.table(opt$normalized_counts, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE)
cat(sprintf("  Loaded matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))

cat("Loading metadata...\n")
metadata <- read.csv(opt$metadata, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Loaded metadata for %d samples\n", nrow(metadata)))

expr_matrix <- as.matrix(counts)

# Check signature gene coverage
cat("\nChecking gene signature coverage in expression data...\n")
available_genes <- rownames(expr_matrix)
for (sig in names(ball_signatures)) {
    overlap <- sum(ball_signatures[[sig]] %in% available_genes)
    total <- length(ball_signatures[[sig]])
    cat(sprintf("  %s: %d/%d genes found (%.1f%%)\n", sig, overlap, total, 100 * overlap / total))
}

# ============================================================================
# Run GSVA/ssGSEA scoring
# ============================================================================
cat("\n--- Running ssGSEA scoring ---\n")

# Filter signatures to only include available genes
filtered_signatures <- lapply(ball_signatures, function(genes) {
    genes[genes %in% available_genes]
})

# Remove signatures with too few genes
min_genes <- 5
valid_sigs <- sapply(filtered_signatures, length) >= min_genes
filtered_signatures <- filtered_signatures[valid_sigs]
cat(sprintf("  %d signatures with >= %d genes available\n",
            sum(valid_sigs), min_genes))

if (length(filtered_signatures) == 0) {
    stop("ERROR: No signatures have sufficient gene coverage. Check gene naming convention.")
}

# Run ssGSEA using GSVA package
cat("  Computing ssGSEA scores...\n")
gsva_param <- ssgseaParam(expr_matrix, filtered_signatures, normalize = TRUE)
ssgsea_scores <- gsva(gsva_param, verbose = FALSE)

cat(sprintf("  ssGSEA scores: %d subtypes x %d samples\n",
            nrow(ssgsea_scores), ncol(ssgsea_scores)))

# Also run GSVA for comparison
cat("  Computing GSVA scores...\n")
gsva_param2 <- gsvaParam(expr_matrix, filtered_signatures, kcdf = "Gaussian")
gsva_scores <- gsva(gsva_param2, verbose = FALSE)

cat(sprintf("  GSVA scores: %d subtypes x %d samples\n",
            nrow(gsva_scores), ncol(gsva_scores)))

# ============================================================================
# Assign subtypes
# ============================================================================
cat("\n--- Assigning molecular subtypes ---\n")

# Use ssGSEA scores for subtype assignment
score_matrix <- as.data.frame(t(ssgsea_scores))
score_matrix$sample_id <- rownames(score_matrix)

# For each sample, find the highest scoring subtype
assignments <- data.frame(
    sample_id = colnames(ssgsea_scores),
    stringsAsFactors = FALSE
)

# Get max score and corresponding subtype for each sample
max_scores <- apply(ssgsea_scores, 2, max)
max_subtypes <- rownames(ssgsea_scores)[apply(ssgsea_scores, 2, which.max)]

# Get second-highest score for confidence margin
second_scores <- apply(ssgsea_scores, 2, function(x) sort(x, decreasing = TRUE)[2])

assignments$assigned_subtype <- max_subtypes
assignments$max_score <- max_scores
assignments$second_score <- second_scores
assignments$score_margin <- max_scores - second_scores

# Apply confidence threshold
assignments$confident <- assignments$max_score >= opt$confidence_threshold
assignments$final_subtype <- ifelse(
    assignments$confident,
    assignments$assigned_subtype,
    "Unclassified"
)

cat(sprintf("  Confident assignments: %d / %d samples (%.1f%%)\n",
            sum(assignments$confident), nrow(assignments),
            100 * sum(assignments$confident) / nrow(assignments)))
cat("\n  Subtype distribution:\n")
print(table(assignments$final_subtype))

# ============================================================================
# Save results
# ============================================================================
cat("\n--- Saving results ---\n")

# Subtype assignments
write.table(assignments %>% select(sample_id, final_subtype, max_score, score_margin, confident),
            file.path(opt$output_dir, "assignments", "subtypes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Full classifier scores (ssGSEA)
classifier_scores <- as.data.frame(t(ssgsea_scores))
classifier_scores$sample_id <- rownames(classifier_scores)
classifier_scores <- classifier_scores %>% select(sample_id, everything())
write.table(classifier_scores,
            file.path(opt$output_dir, "scores", "classifier_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# GSVA scores
gsva_scores_df <- as.data.frame(t(gsva_scores))
gsva_scores_df$sample_id <- rownames(gsva_scores_df)
gsva_scores_df <- gsva_scores_df %>% select(sample_id, everything())
write.table(gsva_scores_df,
            file.path(opt$output_dir, "scores", "gsva_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("  Saved subtypes.tsv, classifier_scores.tsv, gsva_scores.tsv\n")

# ============================================================================
# Generate Plots
# ============================================================================
cat("\n--- Generating Plots ---\n")
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --- Heatmap of ssGSEA scores ---
cat("  Creating score heatmap...\n")
tryCatch({
    # Annotation for assigned subtype
    annotation_col <- data.frame(
        Subtype = assignments$final_subtype,
        row.names = assignments$sample_id
    )

    # Define colors for subtypes
    subtype_colors <- c(
        "ETV6_RUNX1" = "#E41A1C", "BCR_ABL1" = "#377EB8", "KMT2A" = "#4DAF4A",
        "Hyperdiploid" = "#984EA3", "Ph_like" = "#FF7F00", "DUX4" = "#FFFF33",
        "iAMP21" = "#A65628", "PAX5alt" = "#F781BF", "MEF2D" = "#999999",
        "ZNF384" = "#66C2A5", "NUTM1" = "#FC8D62", "Unclassified" = "#CCCCCC"
    )
    # Only use colors for subtypes present
    present_subtypes <- unique(assignments$final_subtype)
    ann_colors <- list(Subtype = subtype_colors[present_subtypes[present_subtypes %in% names(subtype_colors)]])

    # Order samples by subtype assignment
    sample_order <- assignments %>%
        arrange(final_subtype, desc(max_score)) %>%
        pull(sample_id)
    score_mat_ordered <- ssgsea_scores[, sample_order]

    pdf(file.path(plot_dir, "subtype_score_heatmap.pdf"), width = 14, height = 8)
    pheatmap(score_mat_ordered,
             main = "B-ALL Molecular Subtype Scores (ssGSEA)",
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             clustering_distance_rows = "euclidean",
             cluster_cols = FALSE,
             clustering_method = "ward.D2",
             show_colnames = (ncol(score_mat_ordered) <= 50),
             fontsize_col = 6,
             fontsize_row = 10,
             annotation_col = annotation_col[sample_order, , drop = FALSE],
             annotation_colors = ann_colors,
             gaps_col = cumsum(table(assignments$final_subtype[match(sample_order, assignments$sample_id)])))
    dev.off()
    cat("    Heatmap saved.\n")
}, error = function(e) {
    cat(sprintf("    WARNING: Heatmap failed: %s\n", conditionMessage(e)))
})

# --- Barplot of subtype assignments ---
cat("  Creating subtype distribution barplot...\n")
tryCatch({
    subtype_counts <- assignments %>%
        count(final_subtype) %>%
        arrange(desc(n))

    p_bar <- ggplot(subtype_counts, aes(x = reorder(final_subtype, -n), y = n,
                                         fill = final_subtype)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.position = "none"
        ) +
        scale_fill_manual(values = subtype_colors) +
        labs(
            title = "B-ALL Molecular Subtype Distribution",
            x = "Subtype",
            y = "Number of Samples"
        )
    ggsave(file.path(plot_dir, "subtype_distribution.pdf"), p_bar, width = 10, height = 6)
    cat("    Barplot saved.\n")
}, error = function(e) {
    cat(sprintf("    WARNING: Barplot failed: %s\n", conditionMessage(e)))
})

# --- Score confidence plot ---
cat("  Creating confidence score plot...\n")
tryCatch({
    p_conf <- ggplot(assignments, aes(x = reorder(sample_id, -max_score),
                                       y = max_score, fill = final_subtype)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = opt$confidence_threshold, linetype = "dashed",
                   color = "red", linewidth = 0.8) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
            legend.position = "right"
        ) +
        scale_fill_manual(values = subtype_colors) +
        labs(
            title = "Subtype Assignment Confidence Scores",
            subtitle = paste("Dashed line = confidence threshold:", opt$confidence_threshold),
            x = "Sample",
            y = "Maximum ssGSEA Score",
            fill = "Assigned Subtype"
        )
    ggsave(file.path(plot_dir, "subtype_confidence_scores.pdf"),
           p_conf, width = 14, height = 6)
    cat("    Confidence plot saved.\n")
}, error = function(e) {
    cat(sprintf("    WARNING: Confidence plot failed: %s\n", conditionMessage(e)))
})

# --- Score margin plot (max - second) ---
cat("  Creating score margin plot...\n")
tryCatch({
    p_margin <- ggplot(assignments, aes(x = max_score, y = score_margin,
                                         color = final_subtype)) +
        geom_point(size = 2, alpha = 0.7) +
        theme_minimal() +
        scale_color_manual(values = subtype_colors) +
        labs(
            title = "Subtype Classification Confidence",
            x = "Max Signature Score",
            y = "Score Margin (Max - Second)",
            color = "Assigned Subtype"
        )
    ggsave(file.path(plot_dir, "subtype_score_margin.pdf"),
           p_margin, width = 10, height = 7)
    cat("    Score margin plot saved.\n")
}, error = function(e) {
    cat(sprintf("    WARNING: Score margin plot failed: %s\n", conditionMessage(e)))
})

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Molecular Subtyping Complete ===\n")
cat(sprintf("Output directory: %s\n", opt$output_dir))
cat(sprintf("Samples classified: %d / %d (%.1f%% confident)\n",
            sum(assignments$confident), nrow(assignments),
            100 * sum(assignments$confident) / nrow(assignments)))
cat("\nSubtype distribution:\n")
print(table(assignments$final_subtype))
cat("\nKey outputs:\n")
cat(sprintf("  - Subtype assignments: %s\n", file.path(opt$output_dir, "assignments", "subtypes.tsv")))
cat(sprintf("  - Classifier scores: %s\n", file.path(opt$output_dir, "scores", "classifier_scores.tsv")))
cat(sprintf("  - Plots: %s\n", plot_dir))
