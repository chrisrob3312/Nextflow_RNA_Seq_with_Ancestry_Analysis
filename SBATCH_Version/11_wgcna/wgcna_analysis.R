#!/usr/bin/env Rscript
# ============================================================================
# WGCNA: Weighted Gene Co-expression Network Analysis
# ============================================================================
# Builds signed co-expression networks, identifies gene modules, computes
# module eigengenes, and tests correlations with clinical/ancestry traits.
# ============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dynamicTreeCut)
})

# Enable multi-threading for WGCNA
allowWGCNAThreads()

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
  make_option("--norm_counts", type = "character",
              help = "Path to normalized count matrix (genes x samples)"),
  make_option("--metadata", type = "character",
              help = "Path to sample metadata file"),
  make_option("--ancestry_proportions", type = "character", default = "none",
              help = "Path to ancestry proportions file"),
  make_option("--output_dir", type = "character",
              help = "Output directory"),
  make_option("--min_module_size", type = "integer", default = 30,
              help = "Minimum module size for dynamic tree cut"),
  make_option("--merge_height", type = "double", default = 0.25,
              help = "Cut height for merging similar modules"),
  make_option("--network_type", type = "character", default = "signed",
              help = "Network type: signed, unsigned, or signed hybrid"),
  make_option("--top_var_genes", type = "integer", default = 5000,
              help = "Number of top-variance genes to use"),
  make_option("--threads", type = "integer", default = 16,
              help = "Number of threads")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required
for (arg in c("norm_counts", "metadata", "output_dir")) {
  if (is.null(opt[[arg]])) stop(sprintf("Required argument --%s is missing", arg))
}

enableWGCNAThreads(nThreads = opt$threads)

cat("============================================================\n")
cat("WGCNA - Weighted Gene Co-expression Network Analysis\n")
cat("============================================================\n")
cat(sprintf("Normalized counts: %s\n", opt$norm_counts))
cat(sprintf("Metadata:          %s\n", opt$metadata))
cat(sprintf("Ancestry:          %s\n", opt$ancestry_proportions))
cat(sprintf("Output dir:        %s\n", opt$output_dir))
cat(sprintf("Network type:      %s\n", opt$network_type))
cat(sprintf("Min module size:   %d\n", opt$min_module_size))
cat(sprintf("Merge height:      %.3f\n", opt$merge_height))
cat(sprintf("Top var genes:     %d\n", opt$top_var_genes))
cat(sprintf("Threads:           %d\n", opt$threads))
cat("============================================================\n\n")

# ============================================================================
# [1] Load and prepare data
# ============================================================================
cat("[1] Loading data...\n")

# Load normalized counts (genes x samples)
norm_counts <- read.delim(opt$norm_counts, row.names = 1, check.names = FALSE)
cat(sprintf("    Input matrix: %d genes x %d samples\n", nrow(norm_counts), ncol(norm_counts)))

# Load metadata
metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE)
if (!"sample_id" %in% colnames(metadata)) {
  colnames(metadata)[1] <- "sample_id"
}
rownames(metadata) <- metadata$sample_id

# Load ancestry proportions
ancestry_cols <- character(0)
if (opt$ancestry_proportions != "none" && file.exists(opt$ancestry_proportions)) {
  cat("    Loading ancestry proportions...\n")
  ancestry <- read.delim(opt$ancestry_proportions, stringsAsFactors = FALSE)
  if (!"sample_id" %in% colnames(ancestry)) {
    colnames(ancestry)[1] <- "sample_id"
  }
  ancestry_cols <- setdiff(colnames(ancestry), c("sample_id", "predicted_ancestry", "population"))
  ancestry_cols <- ancestry_cols[sapply(ancestry[, ancestry_cols, drop = FALSE], is.numeric)]

  metadata <- merge(metadata, ancestry[, c("sample_id", ancestry_cols)],
                    by = "sample_id", all.x = TRUE)
  rownames(metadata) <- metadata$sample_id
  cat(sprintf("    Ancestry components: %s\n", paste(ancestry_cols, collapse = ", ")))
}

# Match samples
common_samples <- intersect(colnames(norm_counts), rownames(metadata))
if (length(common_samples) < 15) {
  stop("ERROR: WGCNA requires at least 15 samples. Found: ", length(common_samples))
}
norm_counts <- norm_counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]
cat(sprintf("    Matched samples: %d\n", length(common_samples)))

# ============================================================================
# [2] Filter genes by variance
# ============================================================================
cat("\n[2] Filtering to top-variance genes...\n")

gene_vars <- apply(norm_counts, 1, var)
gene_vars <- sort(gene_vars, decreasing = TRUE)

n_genes <- min(opt$top_var_genes, length(gene_vars))
top_genes <- names(gene_vars)[1:n_genes]
expr_data <- norm_counts[top_genes, , drop = FALSE]
cat(sprintf("    Selected top %d genes by variance\n", n_genes))
cat(sprintf("    Variance range: %.3f - %.3f\n", min(gene_vars[top_genes]), max(gene_vars[top_genes])))

# Transpose for WGCNA (samples as rows, genes as columns)
datExpr <- t(expr_data)

# ============================================================================
# [3] Check sample quality
# ============================================================================
cat("\n[3] Checking sample quality...\n")

# Hierarchical clustering to detect outliers
sampleTree <- hclust(dist(datExpr), method = "average")

# Detect outliers using the standardized connectivity method
sample_connectivity <- softConnectivity(datExpr, power = 6, type = opt$network_type)
z_connectivity <- scale(sample_connectivity)
outlier_samples <- rownames(datExpr)[abs(z_connectivity) > 3]

if (length(outlier_samples) > 0) {
  cat(sprintf("    WARNING: %d potential outlier samples detected:\n", length(outlier_samples)))
  cat(sprintf("      %s\n", paste(outlier_samples, collapse = ", ")))
  cat("    These samples are retained but flagged.\n")
} else {
  cat("    No outlier samples detected.\n")
}

# Save sample dendrogram
pdf(file.path(opt$output_dir, "plots", "sample_dendrogram.pdf"), width = 12, height = 6)
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
cat("    Saved: sample_dendrogram.pdf\n")

# ============================================================================
# [4] Determine soft-thresholding power
# ============================================================================
cat("\n[4] Determining soft-thresholding power...\n")

powers <- c(1:20)
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = opt$network_type,
  verbose = 0
)

# Find the lowest power where scale-free topology fit R^2 > 0.85
r2_threshold <- 0.85
sft_df <- sft$fitIndices
valid_powers <- sft_df$Power[sft_df$SFT.R.sq >= r2_threshold]

if (length(valid_powers) > 0) {
  soft_power <- min(valid_powers)
  cat(sprintf("    Selected soft power: %d (R^2 = %.3f)\n",
              soft_power, sft_df$SFT.R.sq[sft_df$Power == soft_power]))
} else {
  # Fall back to the power with highest R^2
  soft_power <- sft_df$Power[which.max(sft_df$SFT.R.sq)]
  cat(sprintf("    WARNING: No power achieved R^2 > %.2f\n", r2_threshold))
  cat(sprintf("    Using power %d with highest R^2 = %.3f\n",
              soft_power, max(sft_df$SFT.R.sq)))
}

# Plot scale-free topology fit
pdf(file.path(opt$output_dir, "plots", "scale_free_topology_fit.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

# Scale-free fit index
plot(sft_df$Power, -sign(sft_df$slope) * sft_df$SFT.R.sq,
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n", main = "Scale Independence")
text(sft_df$Power, -sign(sft_df$slope) * sft_df$SFT.R.sq,
     labels = powers, cex = 0.9, col = "red")
abline(h = r2_threshold, col = "blue", lty = 2)
abline(v = soft_power, col = "darkgreen", lty = 2)

# Mean connectivity
plot(sft_df$Power, sft_df$mean.k.,
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft_df$Power, sft_df$mean.k., labels = powers, cex = 0.9, col = "red")
abline(v = soft_power, col = "darkgreen", lty = 2)
dev.off()
cat("    Saved: scale_free_topology_fit.pdf\n")

# ============================================================================
# [5] Build network and compute TOM
# ============================================================================
cat("\n[5] Building network (this may take a while)...\n")
cat(sprintf("    Power: %d, Network type: %s\n", soft_power, opt$network_type))

# Calculate adjacency matrix
adjacency <- adjacency(datExpr, power = soft_power, type = opt$network_type)

# Calculate Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency, TOMType = opt$network_type)
dissTOM <- 1 - TOM
colnames(dissTOM) <- colnames(datExpr)
rownames(dissTOM) <- colnames(datExpr)

cat("    Adjacency and TOM computed.\n")

# ============================================================================
# [6] Module detection via dynamic tree cut
# ============================================================================
cat("\n[6] Detecting modules...\n")

# Hierarchical clustering of genes using TOM dissimilarity
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Dynamic tree cut for module identification
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = opt$min_module_size
)

# Convert to colors
dynamicColors <- labels2colors(dynamicMods)
cat(sprintf("    Modules before merging: %d\n", length(unique(dynamicColors)) - 1))

# ============================================================================
# [7] Merge similar modules
# ============================================================================
cat("\n[7] Merging similar modules...\n")

# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

# Merge modules with correlation > (1 - merge_height)
merge <- mergeCloseModules(datExpr, dynamicColors,
                            cutHeight = opt$merge_height, verbose = 0)
moduleColors <- merge$colors
MEs <- merge$newMEs

n_modules <- length(unique(moduleColors)) - as.integer("grey" %in% unique(moduleColors))
cat(sprintf("    Modules after merging: %d\n", n_modules))
cat(sprintf("    Module sizes:\n"))
mod_table <- sort(table(moduleColors), decreasing = TRUE)
for (mod in names(mod_table)) {
  if (mod != "grey") {
    cat(sprintf("      %s: %d genes\n", mod, mod_table[mod]))
  }
}
cat(sprintf("      grey (unassigned): %d genes\n",
            sum(moduleColors == "grey")))

# Save module assignments
module_df <- data.frame(
  gene = colnames(datExpr),
  module = moduleColors,
  stringsAsFactors = FALSE
)
write.table(module_df, file = file.path(opt$output_dir, "modules", "module_assignments.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# [8] Plot gene dendrogram with module colors
# ============================================================================
cat("\n[8] Generating dendrogram plot...\n")

pdf(file.path(opt$output_dir, "plots", "gene_dendrogram_modules.pdf"), width = 14, height = 8)
plotDendroAndColors(
  geneTree,
  cbind(dynamicColors, moduleColors),
  c("Dynamic Cut", "Merged Modules"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Dendrogram and Module Colors"
)
dev.off()
cat("    Saved: gene_dendrogram_modules.pdf\n")

# ============================================================================
# [9] Module eigengenes and trait correlations
# ============================================================================
cat("\n[9] Computing module-trait correlations...\n")

# Prepare trait data (numeric traits only)
trait_cols <- character(0)

# Include ancestry proportions
trait_cols <- c(trait_cols, ancestry_cols)

# Include other numeric metadata columns
numeric_cols <- colnames(metadata)[sapply(metadata, is.numeric)]
numeric_cols <- setdiff(numeric_cols, c("sample_id"))
trait_cols <- unique(c(trait_cols, numeric_cols))

# Filter to columns that exist and have variance
trait_data <- metadata[, trait_cols[trait_cols %in% colnames(metadata)], drop = FALSE]
trait_data <- trait_data[, sapply(trait_data, function(x) var(x, na.rm = TRUE) > 0), drop = FALSE]

if (ncol(trait_data) > 0) {
  cat(sprintf("    Trait variables: %s\n", paste(colnames(trait_data), collapse = ", ")))

  # Compute correlations between module eigengenes and traits
  nSamples <- nrow(datExpr)
  moduleTraitCor <- cor(MEs, trait_data, use = "pairwise.complete.obs")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

  # Save correlation results
  cor_df <- as.data.frame(moduleTraitCor)
  cor_df$module <- rownames(cor_df)
  write.table(cor_df, file = file.path(opt$output_dir, "modules", "module_trait_correlations.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  pval_df <- as.data.frame(moduleTraitPvalue)
  pval_df$module <- rownames(pval_df)
  write.table(pval_df, file = file.path(opt$output_dir, "modules", "module_trait_pvalues.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Generate module-trait heatmap
  cat("    Generating module-trait heatmap...\n")

  # Text matrix with correlation and p-value
  textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")")
  dim(textMatrix) <- dim(moduleTraitCor)

  pdf(file.path(opt$output_dir, "plots", "module_trait_heatmap.pdf"),
      width = max(8, ncol(trait_data) * 1.2), height = max(6, nrow(MEs) * 0.5))
  par(mar = c(8, 10, 3, 3))
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = colnames(trait_data),
    yLabels = colnames(MEs),
    ySymbols = gsub("^ME", "", colnames(MEs)),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.6,
    zlim = c(-1, 1),
    main = "Module-Trait Relationships"
  )
  dev.off()
  cat("    Saved: module_trait_heatmap.pdf\n")

  # Also create a pheatmap version (cleaner for many traits)
  pdf(file.path(opt$output_dir, "plots", "module_trait_heatmap_pheatmap.pdf"),
      width = max(8, ncol(trait_data) * 0.8), height = max(6, nrow(MEs) * 0.4))
  # Mark significant correlations
  sig_matrix <- moduleTraitCor
  sig_matrix[moduleTraitPvalue > 0.05] <- NA

  pheatmap(
    moduleTraitCor,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = matrix(ifelse(moduleTraitPvalue < 0.001, "***",
                                     ifelse(moduleTraitPvalue < 0.01, "**",
                                            ifelse(moduleTraitPvalue < 0.05, "*", ""))),
                             nrow = nrow(moduleTraitPvalue)),
    fontsize_number = 8,
    main = "Module-Trait Correlations",
    labels_row = gsub("^ME", "", colnames(MEs))
  )
  dev.off()
  cat("    Saved: module_trait_heatmap_pheatmap.pdf\n")

  # Report significant correlations with ancestry
  if (length(ancestry_cols) > 0) {
    cat("\n    Significant module-ancestry correlations (p < 0.05):\n")
    for (anc in ancestry_cols) {
      if (anc %in% colnames(moduleTraitCor)) {
        sig_mods <- which(moduleTraitPvalue[, anc] < 0.05)
        if (length(sig_mods) > 0) {
          for (idx in sig_mods) {
            cat(sprintf("      %s ~ %s: r=%.3f, p=%.2e\n",
                        gsub("^ME", "", rownames(moduleTraitCor)[idx]),
                        anc,
                        moduleTraitCor[idx, anc],
                        moduleTraitPvalue[idx, anc]))
          }
        }
      }
    }
  }
} else {
  cat("    No numeric traits available for correlation analysis.\n")
}

# ============================================================================
# [10] Identify hub genes per module
# ============================================================================
cat("\n[10] Identifying hub genes...\n")

# Module membership (correlation of gene expression with module eigengene)
geneModuleMembership <- cor(datExpr, MEs, use = "pairwise.complete.obs")
MMPvalue <- corPvalueStudent(geneModuleMembership, nrow(datExpr))

# For each module, identify hub genes (top intramodular connectivity)
hub_gene_list <- list()

for (mod in unique(moduleColors)) {
  if (mod == "grey") next

  mod_genes <- colnames(datExpr)[moduleColors == mod]
  if (length(mod_genes) < 2) next

  # Intramodular connectivity
  mod_idx <- which(moduleColors == mod)
  intra_connectivity <- colSums(adjacency[mod_idx, mod_idx])
  names(intra_connectivity) <- colnames(datExpr)[mod_idx]

  # Module membership for this module
  me_col <- paste0("ME", mod)
  if (me_col %in% colnames(geneModuleMembership)) {
    mm_values <- geneModuleMembership[mod_genes, me_col]
  } else {
    mm_values <- rep(NA, length(mod_genes))
    names(mm_values) <- mod_genes
  }

  # Combine metrics
  hub_df <- data.frame(
    gene = mod_genes,
    module = mod,
    intramodular_connectivity = intra_connectivity[mod_genes],
    module_membership = mm_values[mod_genes],
    stringsAsFactors = FALSE
  )
  hub_df <- hub_df[order(-hub_df$intramodular_connectivity), ]
  hub_df$hub_rank <- seq_len(nrow(hub_df))

  # Save per-module hub genes
  write.table(hub_df,
              file = file.path(opt$output_dir, "hub_genes", sprintf("hub_genes_%s.tsv", mod)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Store top hubs
  hub_gene_list[[mod]] <- head(hub_df, 20)
}

# Save combined hub genes summary (top 10 per module)
all_hubs <- do.call(rbind, lapply(hub_gene_list, function(x) head(x, 10)))
write.table(all_hubs, file = file.path(opt$output_dir, "hub_genes", "all_hub_genes_top10.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("    Hub genes identified for %d modules\n", length(hub_gene_list)))

# ============================================================================
# [11] Save eigengenes and connectivity
# ============================================================================
cat("\n[11] Saving module eigengenes and network properties...\n")

# Save module eigengenes
ME_df <- as.data.frame(MEs)
ME_df$sample_id <- rownames(ME_df)
ME_df <- ME_df[, c("sample_id", setdiff(colnames(ME_df), "sample_id"))]
write.table(ME_df, file = file.path(opt$output_dir, "modules", "module_eigengenes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Save gene-level network properties
gene_props <- data.frame(
  gene = colnames(datExpr),
  module = moduleColors,
  total_connectivity = colSums(adjacency),
  stringsAsFactors = FALSE
)

# Add module membership for assigned module
gene_props$module_membership <- sapply(seq_len(nrow(gene_props)), function(i) {
  me_col <- paste0("ME", gene_props$module[i])
  if (me_col %in% colnames(geneModuleMembership)) {
    return(geneModuleMembership[gene_props$gene[i], me_col])
  }
  return(NA)
})

write.table(gene_props, file = file.path(opt$output_dir, "network", "gene_network_properties.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# [12] Save R objects for downstream use
# ============================================================================
cat("\n[12] Saving R objects...\n")

save(datExpr, moduleColors, MEs, geneTree, TOM, adjacency, soft_power,
     geneModuleMembership, MMPvalue, module_df,
     file = file.path(opt$output_dir, "rdata", "wgcna_network.RData"))
cat("    Saved: wgcna_network.RData\n")

# Save soft power info
save(sft, soft_power,
     file = file.path(opt$output_dir, "rdata", "wgcna_soft_power.RData"))

cat("\n============================================================\n")
cat("WGCNA analysis complete!\n")
cat(sprintf("  Modules detected: %d\n", n_modules))
cat(sprintf("  Soft power used: %d\n", soft_power))
cat(sprintf("  Total genes in network: %d\n", ncol(datExpr)))
cat(sprintf("  Unassigned (grey) genes: %d\n", sum(moduleColors == "grey")))
cat("============================================================\n")

# Session info
writeLines(capture.output(sessionInfo()),
           file.path(opt$output_dir, "session_info.txt"))
