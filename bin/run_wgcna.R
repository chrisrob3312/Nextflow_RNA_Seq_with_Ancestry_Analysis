#!/usr/bin/env Rscript

# ============================================================================
# WGCNA - Weighted Gene Co-expression Network Analysis
# Correlates modules with traits: ancestry, relapse, ADI, cytogenetics
# ============================================================================

suppressPackageStartupMessages({
    library(WGCNA)
    library(optparse)
    library(ggplot2)
    library(pheatmap)
    library(dynamicTreeCut)
})

allowWGCNAThreads()

option_list <- list(
    make_option("--expression", type = "character", help = "Normalized expression matrix"),
    make_option("--metadata", type = "character", help = "Sample metadata TSV"),
    make_option("--ancestry", type = "character", help = "Ancestry proportions TSV"),
    make_option("--min-module-size", type = "integer", default = 30),
    make_option("--merge-cut-height", type = "double", default = 0.25),
    make_option("--network-type", type = "character", default = "signed"),
    make_option("--tom-type", type = "character", default = "signed"),
    make_option("--soft-power", type = "integer", default = NULL),
    make_option("--output-dir", type = "character", default = "wgcna_results"),
    make_option("--plot-dir", type = "character", default = "wgcna_plots"),
    make_option("--threads", type = "integer", default = 4)
)
opt <- parse_args(OptionParser(option_list = option_list))

enableWGCNAThreads(nThreads = opt$threads)

# ---- Load data ----
expr <- read.delim(opt$expression, row.names = 1, check.names = FALSE)
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

if (!is.null(opt$ancestry) && file.exists(opt$ancestry)) {
    ancestry <- read.delim(opt$ancestry, check.names = FALSE)
    rownames(ancestry) <- ancestry$sample_id
    common <- intersect(rownames(metadata), rownames(ancestry))
    ancestry_cols <- setdiff(colnames(ancestry), "sample_id")
    metadata[common, ancestry_cols] <- ancestry[common, ancestry_cols]
}

common <- intersect(colnames(expr), rownames(metadata))
expr <- expr[, common]
metadata <- metadata[common, ]

# Transpose for WGCNA (samples x genes)
datExpr <- t(expr)

# ---- Filter genes with high variance ----
gene_vars <- apply(datExpr, 2, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(5000, length(gene_vars))]
datExpr <- datExpr[, top_genes]

cat("Using", ncol(datExpr), "genes for WGCNA\n")

# ---- Check for outlier samples ----
sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file.path(opt$`plot-dir`, "sample_dendrogram.pdf"), width = 12, height = 6)
plot(sampleTree, main = "Sample Dendrogram", sub = "", xlab = "")
dev.off()

# ---- Pick soft-thresholding power ----
if (is.null(opt$`soft-power`)) {
    powers <- c(1:20)
    sft <- pickSoftThreshold(datExpr, powerVector = powers,
                              networkType = opt$`network-type`, verbose = 0)

    # Select power where scale-free R^2 > 0.85
    power <- sft$powerEstimate
    if (is.na(power)) power <- 6  # Fallback

    pdf(file.path(opt$`plot-dir`, "soft_threshold.pdf"), width = 10, height = 5)
    par(mfrow = c(1, 2))
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
         xlab = "Soft Threshold (power)", ylab = "Scale Free R^2",
         main = "Scale Independence")
    abline(h = 0.85, col = "red")
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
         xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
         main = "Mean Connectivity")
    dev.off()

    cat("Selected soft-thresholding power:", power, "\n")
} else {
    power <- opt$`soft-power`
}

# ---- Build network and detect modules ----
net <- blockwiseModules(
    datExpr,
    power = power,
    networkType = opt$`network-type`,
    TOMType = opt$`tom-type`,
    minModuleSize = opt$`min-module-size`,
    reassignThreshold = 0,
    mergeCutHeight = opt$`merge-cut-height`,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 3,
    maxBlockSize = ncol(datExpr)
)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

cat("Detected", length(unique(moduleColors)) - 1, "modules (excluding grey)\n")

# Module dendrogram
pdf(file.path(opt$`plot-dir`, "module_dendrogram.pdf"), width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# ---- Module eigengenes ----
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs)

write.table(cbind(sample_id = rownames(MEs), MEs),
            file.path(opt$`output-dir`, "module_eigengenes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Module-trait correlations ----
# Build numeric trait matrix
trait_cols <- c("age", "sex", "blast_percentage", "tumor_purity",
                "relapse_status", "adi_quartile", "timepoint")
ancestry_cols <- grep("^pct_", colnames(metadata), value = TRUE)
trait_cols <- c(trait_cols, ancestry_cols)
trait_cols <- trait_cols[trait_cols %in% colnames(metadata)]

trait_data <- metadata[rownames(MEs), trait_cols, drop = FALSE]

# Encode categorical variables numerically
for (col in colnames(trait_data)) {
    if (is.character(trait_data[[col]]) || is.factor(trait_data[[col]])) {
        vals <- trait_data[[col]]
        if (all(vals %in% c("M", "F", "NA"))) {
            trait_data[[col]] <- ifelse(vals == "M", 1, ifelse(vals == "F", 0, NA))
        } else if (all(vals %in% c("relapse", "no_relapse", "NA"))) {
            trait_data[[col]] <- ifelse(vals == "relapse", 1, ifelse(vals == "no_relapse", 0, NA))
        } else if (all(vals %in% c("diagnostic", "relapse", "NA"))) {
            trait_data[[col]] <- ifelse(vals == "relapse", 1, ifelse(vals == "diagnostic", 0, NA))
        } else {
            trait_data[[col]] <- as.numeric(as.factor(trait_data[[col]]))
        }
    }
    trait_data[[col]] <- as.numeric(trait_data[[col]])
}

# Correlation
module_trait_cor <- cor(MEs, trait_data, use = "pairwise.complete.obs")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

write.table(module_trait_cor, file.path(opt$`output-dir`, "module_trait_cor.tsv"),
            sep = "\t", quote = FALSE)

# Heatmap
pdf(file.path(opt$`plot-dir`, "module_trait_heatmap.pdf"), width = 12, height = 10)
textMatrix <- paste0(signif(module_trait_cor, 2), "\n(",
                      signif(module_trait_pval, 1), ")")
dim(textMatrix) <- dim(module_trait_cor)
par(mar = c(8, 10, 3, 3))
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = colnames(trait_data),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-Trait Correlations")
dev.off()

# ---- Module membership and hub genes ----
geneModuleMembership <- cor(datExpr, MEs, use = "p")
module_membership_df <- data.frame(
    gene = colnames(datExpr),
    module = moduleColors,
    stringsAsFactors = FALSE
)
for (me_col in colnames(MEs)) {
    module_membership_df[[paste0("MM_", me_col)]] <- geneModuleMembership[, me_col]
}

write.table(module_membership_df,
            file.path(opt$`output-dir`, "module_membership.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Hub genes (top connectivity per module)
hub_genes <- data.frame()
for (mod in unique(moduleColors[moduleColors != "grey"])) {
    mod_genes <- colnames(datExpr)[moduleColors == mod]
    me_col <- paste0("ME", mod)
    if (me_col %in% colnames(MEs)) {
        mm <- abs(geneModuleMembership[mod_genes, me_col])
        top_hubs <- head(sort(mm, decreasing = TRUE), 20)
        hub_genes <- rbind(hub_genes, data.frame(
            module = mod, gene = names(top_hubs),
            module_membership = as.numeric(top_hubs)
        ))
    }
}
write.table(hub_genes, file.path(opt$`output-dir`, "hub_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("WGCNA analysis complete.\n")
