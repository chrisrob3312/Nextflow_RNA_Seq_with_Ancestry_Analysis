#!/usr/bin/env Rscript
# ============================================================================
# Differential Expression Analysis: DESeq2 + limma-voom
# ============================================================================
# Runs both DESeq2 and limma-voom on all contrasts defined in a JSON file.
# Supports categorical, continuous, and categorical_multi contrast types.
# Incorporates ancestry proportions and user-defined covariates.
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(jsonlite)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(RColorBrewer)
  library(BiocParallel)
  library(optparse)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
  make_option("--count_matrix", type = "character", help = "Path to raw count matrix (genes x samples)"),
  make_option("--metadata", type = "character", help = "Path to sample metadata file"),
  make_option("--ancestry_proportions", type = "character", default = "none",
              help = "Path to ancestry proportions file (or 'none')"),
  make_option("--contrasts_json", type = "character", help = "Path to contrasts JSON file"),
  make_option("--covariates", type = "character", default = "",
              help = "Comma-separated covariates to include in the model"),
  make_option("--output_dir", type = "character", help = "Output directory"),
  make_option("--padj_threshold", type = "double", default = 0.05,
              help = "Adjusted p-value threshold for significance"),
  make_option("--lfc_threshold", type = "double", default = 0.585,
              help = "Log2 fold change threshold for significance"),
  make_option("--threads", type = "integer", default = 8,
              help = "Number of threads for parallel processing")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
required_args <- c("count_matrix", "metadata", "contrasts_json", "output_dir")
for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(sprintf("Required argument --%s is missing", arg))
  }
}

# Set up parallel processing
register(MulticoreParam(workers = opt$threads))

cat("============================================================\n")
cat("Differential Expression Analysis\n")
cat("============================================================\n")
cat(sprintf("Count matrix:   %s\n", opt$count_matrix))
cat(sprintf("Metadata:       %s\n", opt$metadata))
cat(sprintf("Ancestry:       %s\n", opt$ancestry_proportions))
cat(sprintf("Contrasts:      %s\n", opt$contrasts_json))
cat(sprintf("Covariates:     %s\n", opt$covariates))
cat(sprintf("Output dir:     %s\n", opt$output_dir))
cat(sprintf("padj cutoff:    %s\n", opt$padj_threshold))
cat(sprintf("LFC cutoff:     %s\n", opt$lfc_threshold))
cat(sprintf("Threads:        %d\n", opt$threads))
cat("============================================================\n\n")

# ============================================================================
# Load and prepare data
# ============================================================================
cat("[1] Loading data...\n")

# Load count matrix
counts <- read.delim(opt$count_matrix, row.names = 1, check.names = FALSE)
cat(sprintf("    Count matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))

# Load metadata
metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE)
if (!"sample_id" %in% colnames(metadata)) {
  # Assume first column is sample ID

  colnames(metadata)[1] <- "sample_id"
}
rownames(metadata) <- metadata$sample_id

# Load and merge ancestry proportions
if (opt$ancestry_proportions != "none" && file.exists(opt$ancestry_proportions)) {
  cat("    Loading ancestry proportions...\n")
  ancestry <- read.delim(opt$ancestry_proportions, stringsAsFactors = FALSE)
  if (!"sample_id" %in% colnames(ancestry)) {
    colnames(ancestry)[1] <- "sample_id"
  }
  # Identify ancestry proportion columns (exclude sample_id and any classification columns)
  ancestry_cols <- setdiff(colnames(ancestry), c("sample_id", "predicted_ancestry", "population"))
  # Keep only numeric ancestry columns
  ancestry_cols <- ancestry_cols[sapply(ancestry[, ancestry_cols, drop = FALSE], is.numeric)]
  cat(sprintf("    Ancestry components: %s\n", paste(ancestry_cols, collapse = ", ")))

  # Merge with metadata
  metadata <- merge(metadata, ancestry[, c("sample_id", ancestry_cols)],
                    by = "sample_id", all.x = TRUE)
  rownames(metadata) <- metadata$sample_id
} else {
  ancestry_cols <- character(0)
  cat("    No ancestry proportions loaded.\n")
}

# Ensure samples match between counts and metadata
common_samples <- intersect(colnames(counts), rownames(metadata))
if (length(common_samples) == 0) {
  stop("ERROR: No overlapping samples between count matrix and metadata")
}
cat(sprintf("    Common samples: %d\n", length(common_samples)))

counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# ============================================================================
# Parse contrasts
# ============================================================================
cat("\n[2] Parsing contrasts...\n")
contrasts_list <- fromJSON(opt$contrasts_json, simplifyVector = FALSE)

# If the JSON has a top-level "contrasts" field, use it
if ("contrasts" %in% names(contrasts_list)) {
  contrasts_list <- contrasts_list$contrasts
}

cat(sprintf("    Number of contrasts: %d\n", length(contrasts_list)))

# ============================================================================
# Parse covariates
# ============================================================================
covariate_vars <- character(0)
if (nchar(opt$covariates) > 0 && opt$covariates != "") {
  covariate_vars <- trimws(unlist(strsplit(opt$covariates, ",")))
  # Filter to covariates that exist in metadata
  available_covs <- covariate_vars[covariate_vars %in% colnames(metadata)]
  missing_covs <- setdiff(covariate_vars, colnames(metadata))
  if (length(missing_covs) > 0) {
    cat(sprintf("    WARNING: Covariates not found in metadata: %s\n",
                paste(missing_covs, collapse = ", ")))
  }
  covariate_vars <- available_covs
}

# Add ancestry proportions as covariates (excluding one to avoid collinearity)
if (length(ancestry_cols) > 1) {
  ancestry_covs <- ancestry_cols[1:(length(ancestry_cols) - 1)]
  covariate_vars <- unique(c(covariate_vars, ancestry_covs))
}

cat(sprintf("    Model covariates: %s\n",
            ifelse(length(covariate_vars) > 0, paste(covariate_vars, collapse = ", "), "none")))

# ============================================================================
# Helper Functions
# ============================================================================

build_design_formula <- function(contrast, covariates, contrast_type) {
  # Build the design formula based on contrast type and covariates
  terms <- character(0)

  # Add covariates
  if (length(covariates) > 0) {
    terms <- c(terms, covariates)
  }

  # Add the contrast variable
  if (contrast_type == "categorical" || contrast_type == "categorical_multi") {
    terms <- c(terms, contrast$variable)
  } else if (contrast_type == "continuous") {
    terms <- c(terms, contrast$variable)
  }

  formula_str <- paste("~", paste(terms, collapse = " + "))
  return(as.formula(formula_str))
}

filter_low_counts <- function(counts_mat, min_count = 10, min_samples = 3) {
  keep <- rowSums(counts_mat >= min_count) >= min_samples
  return(counts_mat[keep, , drop = FALSE])
}

generate_volcano_plot <- function(results_df, contrast_name, method, output_dir,
                                  padj_thresh, lfc_thresh) {
  # Create volcano plot using EnhancedVolcano
  p <- EnhancedVolcano(
    results_df,
    lab = rownames(results_df),
    x = "log2FoldChange",
    y = "padj",
    pCutoff = padj_thresh,
    FCcutoff = lfc_thresh,
    title = sprintf("%s - %s", contrast_name, method),
    subtitle = sprintf("padj < %s, |LFC| > %s", padj_thresh, lfc_thresh),
    legendPosition = "right",
    labSize = 3,
    drawConnectors = TRUE,
    widthConnectors = 0.5
  )

  filename <- file.path(output_dir, "plots",
                        sprintf("volcano_%s_%s.pdf", method, gsub("[^A-Za-z0-9_]", "_", contrast_name)))
  ggsave(filename, p, width = 10, height = 8)
  cat(sprintf("      Saved: %s\n", basename(filename)))
}

generate_ma_plot <- function(results_df, contrast_name, method, output_dir,
                             padj_thresh, lfc_thresh) {
  df <- data.frame(
    baseMean = results_df$baseMean,
    log2FoldChange = results_df$log2FoldChange,
    significant = ifelse(!is.na(results_df$padj) &
                           results_df$padj < padj_thresh &
                           abs(results_df$log2FoldChange) > lfc_thresh,
                         "Significant", "NS"),
    stringsAsFactors = FALSE
  )

  p <- ggplot(df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_color_manual(values = c("NS" = "grey60", "Significant" = "firebrick")) +
    geom_hline(yintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "blue") +
    labs(title = sprintf("MA Plot: %s (%s)", contrast_name, method),
         x = "log10(baseMean + 1)", y = "log2 Fold Change") +
    theme_bw() +
    theme(legend.position = "bottom")

  filename <- file.path(output_dir, "plots",
                        sprintf("ma_%s_%s.pdf", method, gsub("[^A-Za-z0-9_]", "_", contrast_name)))
  ggsave(filename, p, width = 8, height = 6)
  cat(sprintf("      Saved: %s\n", basename(filename)))
}

generate_heatmap <- function(norm_counts, results_df, metadata_sub, contrast_name,
                             method, output_dir, n_top = 50) {
  # Select top significant genes by adjusted p-value
  sig_genes <- results_df[order(results_df$padj), , drop = FALSE]
  sig_genes <- sig_genes[!is.na(sig_genes$padj), , drop = FALSE]
  top_genes <- head(rownames(sig_genes), n_top)

  if (length(top_genes) < 2) {
    cat("      Skipping heatmap: fewer than 2 significant genes\n")
    return(invisible(NULL))
  }

  # Subset and scale
  mat <- norm_counts[top_genes, , drop = FALSE]
  mat_scaled <- t(scale(t(mat)))

  # Annotation
  anno_cols <- intersect(c("condition", "group", "predicted_ancestry"), colnames(metadata_sub))
  if (length(anno_cols) > 0) {
    annotation_col <- metadata_sub[, anno_cols, drop = FALSE]
  } else {
    annotation_col <- NA
  }

  filename <- file.path(output_dir, "plots",
                        sprintf("heatmap_%s_%s.pdf", method, gsub("[^A-Za-z0-9_]", "_", contrast_name)))

  pdf(filename, width = 12, height = 10)
  pheatmap(
    mat_scaled,
    annotation_col = annotation_col,
    show_rownames = (length(top_genes) <= 50),
    show_colnames = FALSE,
    clustering_method = "ward.D2",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = sprintf("Top %d DE Genes: %s (%s)", length(top_genes), contrast_name, method),
    fontsize_row = 6
  )
  dev.off()
  cat(sprintf("      Saved: %s\n", basename(filename)))
}

# ============================================================================
# Run DESeq2 analysis for a single contrast
# ============================================================================
run_deseq2 <- function(counts_mat, metadata_sub, contrast, covariates, output_dir,
                       padj_thresh, lfc_thresh) {

  contrast_name <- contrast$name
  contrast_type <- contrast$type
  cat(sprintf("    [DESeq2] Running contrast: %s (type: %s)\n", contrast_name, contrast_type))

  # Build formula
  design_formula <- build_design_formula(contrast, covariates, contrast_type)
  cat(sprintf("      Design: %s\n", deparse(design_formula)))

  # Prepare metadata - ensure factor levels for categorical

  if (contrast_type %in% c("categorical", "categorical_multi")) {
    metadata_sub[[contrast$variable]] <- factor(metadata_sub[[contrast$variable]])
    if (!is.null(contrast$reference)) {
      metadata_sub[[contrast$variable]] <- relevel(metadata_sub[[contrast$variable]],
                                                    ref = contrast$reference)
    }
  }

  # Ensure covariate columns are appropriate types
  for (cov in covariates) {
    if (cov %in% colnames(metadata_sub)) {
      if (is.character(metadata_sub[[cov]])) {
        metadata_sub[[cov]] <- factor(metadata_sub[[cov]])
      }
    }
  }

  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = metadata_sub,
    design = design_formula
  )

  # Filter low-count genes
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  cat(sprintf("      Genes after filtering: %d\n", nrow(dds)))

  # Run DESeq2
  dds <- DESeq(dds, parallel = TRUE)

  # Extract results based on contrast type
  if (contrast_type == "categorical") {
    res <- results(dds,
                   contrast = c(contrast$variable, contrast$target, contrast$reference),
                   alpha = padj_thresh)
  } else if (contrast_type == "categorical_multi") {
    # For multi-level: extract results for each comparison vs reference
    levels_to_test <- setdiff(levels(metadata_sub[[contrast$variable]]), contrast$reference)
    res_list <- list()
    for (lvl in levels_to_test) {
      res_list[[lvl]] <- results(dds,
                                  contrast = c(contrast$variable, lvl, contrast$reference),
                                  alpha = padj_thresh)
    }
    # Return first comparison for main output, save all
    res <- res_list[[1]]
    for (lvl in names(res_list)) {
      sub_name <- sprintf("%s_%s_vs_%s", contrast_name, lvl, contrast$reference)
      res_df <- as.data.frame(res_list[[lvl]])
      res_df <- res_df[order(res_df$padj), ]
      write.table(res_df,
                  file = file.path(output_dir, "deseq2",
                                   sprintf("%s.tsv", gsub("[^A-Za-z0-9_]", "_", sub_name))),
                  sep = "\t", quote = FALSE, col.names = NA)
    }
  } else if (contrast_type == "continuous") {
    res <- results(dds, name = contrast$variable, alpha = padj_thresh)
  }

  # Convert to data frame
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$padj), ]

  # Add gene names as column
  res_df$gene <- rownames(res_df)

  # Summary statistics
  n_sig <- sum(!is.na(res_df$padj) & res_df$padj < padj_thresh &
                 abs(res_df$log2FoldChange) > lfc_thresh, na.rm = TRUE)
  n_up <- sum(!is.na(res_df$padj) & res_df$padj < padj_thresh &
                res_df$log2FoldChange > lfc_thresh, na.rm = TRUE)
  n_down <- sum(!is.na(res_df$padj) & res_df$padj < padj_thresh &
                  res_df$log2FoldChange < -lfc_thresh, na.rm = TRUE)
  cat(sprintf("      Significant: %d (up: %d, down: %d)\n", n_sig, n_up, n_down))

  # Save results table
  out_file <- file.path(output_dir, "deseq2",
                        sprintf("%s.tsv", gsub("[^A-Za-z0-9_]", "_", contrast_name)))
  write.table(res_df, file = out_file, sep = "\t", quote = FALSE, col.names = NA)
  cat(sprintf("      Results saved: %s\n", basename(out_file)))

  # Generate plots
  generate_volcano_plot(res_df, contrast_name, "deseq2", output_dir, padj_thresh, lfc_thresh)
  generate_ma_plot(res_df, contrast_name, "deseq2", output_dir, padj_thresh, lfc_thresh)

  # Heatmap with normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  generate_heatmap(norm_counts, res_df, metadata_sub, contrast_name, "deseq2",
                   output_dir, n_top = 50)

  # Return DESeq2 object and results for downstream use
  return(list(dds = dds, results = res_df, normalized_counts = norm_counts))
}

# ============================================================================
# Run limma-voom analysis for a single contrast
# ============================================================================
run_limma_voom <- function(counts_mat, metadata_sub, contrast, covariates, output_dir,
                           padj_thresh, lfc_thresh) {

  contrast_name <- contrast$name
  contrast_type <- contrast$type
  cat(sprintf("    [limma-voom] Running contrast: %s (type: %s)\n", contrast_name, contrast_type))

  # Prepare metadata - ensure factor levels for categorical

  if (contrast_type %in% c("categorical", "categorical_multi")) {
    metadata_sub[[contrast$variable]] <- factor(metadata_sub[[contrast$variable]])
    if (!is.null(contrast$reference)) {
      metadata_sub[[contrast$variable]] <- relevel(metadata_sub[[contrast$variable]],
                                                    ref = contrast$reference)
    }
  }

  # Ensure covariate columns are appropriate types
  for (cov in covariates) {
    if (cov %in% colnames(metadata_sub)) {
      if (is.character(metadata_sub[[cov]])) {
        metadata_sub[[cov]] <- factor(metadata_sub[[cov]])
      }
    }
  }

  # Build design matrix
  design_formula <- build_design_formula(contrast, covariates, contrast_type)
  design <- model.matrix(design_formula, data = metadata_sub)
  cat(sprintf("      Design matrix: %d samples x %d coefficients\n", nrow(design), ncol(design)))

  # Create DGEList and filter

  dge <- DGEList(counts = counts_mat)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat(sprintf("      Genes after filtering: %d\n", nrow(dge)))

  # TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")

  # Voom transformation with precision weights
  v <- voom(dge, design, plot = FALSE)

  # Fit linear model

  fit <- lmFit(v, design)

  # Set up contrasts based on type

  if (contrast_type == "categorical") {
    coef_name <- sprintf("%s%s", contrast$variable, contrast$target)
    # Find the coefficient matching the contrast
    matching_coefs <- grep(contrast$variable, colnames(design), value = TRUE)
    target_coef <- paste0(contrast$variable, contrast$target)
    if (target_coef %in% colnames(design)) {
      coef_idx <- which(colnames(design) == target_coef)
    } else {
      # Try to find it with a different naming
      coef_idx <- which(grepl(contrast$target, colnames(design)))
      if (length(coef_idx) == 0) coef_idx <- ncol(design)
    }
    fit2 <- eBayes(fit)
    res <- topTable(fit2, coef = coef_idx, number = Inf, sort.by = "none")

  } else if (contrast_type == "categorical_multi") {
    fit2 <- eBayes(fit)
    levels_to_test <- setdiff(levels(metadata_sub[[contrast$variable]]), contrast$reference)

    for (lvl in levels_to_test) {
      target_coef <- paste0(contrast$variable, lvl)
      if (target_coef %in% colnames(design)) {
        coef_idx <- which(colnames(design) == target_coef)
        sub_res <- topTable(fit2, coef = coef_idx, number = Inf, sort.by = "none")
        sub_name <- sprintf("%s_%s_vs_%s", contrast_name, lvl, contrast$reference)

        # Standardize column names
        sub_res_out <- data.frame(
          gene = rownames(sub_res),
          log2FoldChange = sub_res$logFC,
          AveExpr = sub_res$AveExpr,
          t_statistic = sub_res$t,
          pvalue = sub_res$P.Value,
          padj = sub_res$adj.P.Val,
          B = sub_res$B,
          row.names = rownames(sub_res)
        )
        write.table(sub_res_out,
                    file = file.path(output_dir, "limma_voom",
                                     sprintf("%s.tsv", gsub("[^A-Za-z0-9_]", "_", sub_name))),
                    sep = "\t", quote = FALSE, col.names = NA)
      }
    }
    # Use first level for main result
    target_coef <- paste0(contrast$variable, levels_to_test[1])
    coef_idx <- which(colnames(design) == target_coef)
    if (length(coef_idx) == 0) coef_idx <- 2
    res <- topTable(fit2, coef = coef_idx, number = Inf, sort.by = "none")

  } else if (contrast_type == "continuous") {
    coef_idx <- which(colnames(design) == contrast$variable)
    if (length(coef_idx) == 0) coef_idx <- ncol(design)
    fit2 <- eBayes(fit)
    res <- topTable(fit2, coef = coef_idx, number = Inf, sort.by = "none")
  }

  # Standardize output columns to match DESeq2 output format
  res_df <- data.frame(
    gene = rownames(res),
    log2FoldChange = res$logFC,
    baseMean = res$AveExpr,
    t_statistic = res$t,
    pvalue = res$P.Value,
    padj = res$adj.P.Val,
    B = res$B,
    row.names = rownames(res),
    stringsAsFactors = FALSE
  )
  res_df <- res_df[order(res_df$padj), ]

  # Summary
  n_sig <- sum(!is.na(res_df$padj) & res_df$padj < padj_thresh &
                 abs(res_df$log2FoldChange) > lfc_thresh, na.rm = TRUE)
  n_up <- sum(!is.na(res_df$padj) & res_df$padj < padj_thresh &
                res_df$log2FoldChange > lfc_thresh, na.rm = TRUE)
  n_down <- sum(!is.na(res_df$padj) & res_df$padj < padj_thresh &
                  res_df$log2FoldChange < -lfc_thresh, na.rm = TRUE)
  cat(sprintf("      Significant: %d (up: %d, down: %d)\n", n_sig, n_up, n_down))

  # Save results
  out_file <- file.path(output_dir, "limma_voom",
                        sprintf("%s.tsv", gsub("[^A-Za-z0-9_]", "_", contrast_name)))
  write.table(res_df, file = out_file, sep = "\t", quote = FALSE, col.names = NA)
  cat(sprintf("      Results saved: %s\n", basename(out_file)))

  # Generate plots
  generate_volcano_plot(res_df, contrast_name, "limma_voom", output_dir, padj_thresh, lfc_thresh)
  generate_ma_plot(res_df, contrast_name, "limma_voom", output_dir, padj_thresh, lfc_thresh)

  # Heatmap using voom-normalized log-CPM
  norm_counts <- v$E
  generate_heatmap(norm_counts, res_df, metadata_sub, contrast_name, "limma_voom",
                   output_dir, n_top = 50)

  return(list(fit = fit2, voom = v, results = res_df))
}

# ============================================================================
# Main Analysis Loop
# ============================================================================
cat("\n[3] Running differential expression analyses...\n")
cat("============================================================\n\n")

# Store all results for downstream use
all_results <- list()

for (i in seq_along(contrasts_list)) {
  contrast <- contrasts_list[[i]]
  contrast_name <- contrast$name
  contrast_type <- ifelse(is.null(contrast$type), "categorical", contrast$type)
  contrast$type <- contrast_type

  cat(sprintf("\n--- Contrast %d/%d: %s ---\n", i, length(contrasts_list), contrast_name))

  # Subset samples if subgroup filtering is specified
  if (!is.null(contrast$subgroup_variable) && !is.null(contrast$subgroup_value)) {
    cat(sprintf("    Filtering to subgroup: %s == %s\n",
                contrast$subgroup_variable, contrast$subgroup_value))
    sub_idx <- metadata[[contrast$subgroup_variable]] == contrast$subgroup_value
    metadata_sub <- metadata[sub_idx, , drop = FALSE]
    counts_sub <- counts[, rownames(metadata_sub), drop = FALSE]
  } else {
    metadata_sub <- metadata
    counts_sub <- counts
  }

  # Remove samples with NA in the contrast variable
  valid_idx <- !is.na(metadata_sub[[contrast$variable]])
  metadata_sub <- metadata_sub[valid_idx, , drop = FALSE]
  counts_sub <- counts_sub[, rownames(metadata_sub), drop = FALSE]

  cat(sprintf("    Samples in analysis: %d\n", nrow(metadata_sub)))

  if (nrow(metadata_sub) < 4) {
    cat("    WARNING: Fewer than 4 samples, skipping this contrast.\n")
    next
  }

  # Filter covariates to those without too many NAs or singular levels
  valid_covariates <- character(0)
  for (cov in covariate_vars) {
    if (cov %in% colnames(metadata_sub)) {
      col_data <- metadata_sub[[cov]]
      n_na <- sum(is.na(col_data))
      if (n_na / length(col_data) > 0.5) {
        cat(sprintf("    Dropping covariate '%s': >50%% missing\n", cov))
        next
      }
      if (is.character(col_data) || is.factor(col_data)) {
        n_levels <- length(unique(na.omit(col_data)))
        if (n_levels < 2) {
          cat(sprintf("    Dropping covariate '%s': only 1 level\n", cov))
          next
        }
      }
      valid_covariates <- c(valid_covariates, cov)
    }
  }

  # Remove samples with NA in covariates
  if (length(valid_covariates) > 0) {
    complete_idx <- complete.cases(metadata_sub[, valid_covariates, drop = FALSE])
    if (sum(complete_idx) < nrow(metadata_sub)) {
      cat(sprintf("    Removing %d samples with missing covariates\n",
                  sum(!complete_idx)))
      metadata_sub <- metadata_sub[complete_idx, , drop = FALSE]
      counts_sub <- counts_sub[, rownames(metadata_sub), drop = FALSE]
    }
  }

  # Run DESeq2
  tryCatch({
    deseq2_result <- run_deseq2(counts_sub, metadata_sub, contrast, valid_covariates,
                                opt$output_dir, opt$padj_threshold, opt$lfc_threshold)
    all_results[[contrast_name]]$deseq2 <- deseq2_result
  }, error = function(e) {
    cat(sprintf("    ERROR in DESeq2 for '%s': %s\n", contrast_name, conditionMessage(e)))
  })

  # Run limma-voom
  tryCatch({
    limma_result <- run_limma_voom(counts_sub, metadata_sub, contrast, valid_covariates,
                                   opt$output_dir, opt$padj_threshold, opt$lfc_threshold)
    all_results[[contrast_name]]$limma_voom <- limma_result
  }, error = function(e) {
    cat(sprintf("    ERROR in limma-voom for '%s': %s\n", contrast_name, conditionMessage(e)))
  })
}

# ============================================================================
# Save R objects for downstream analysis
# ============================================================================
cat("\n[4] Saving R objects...\n")

save(all_results,
     file = file.path(opt$output_dir, "rdata", "de_all_results.RData"))
cat(sprintf("    Saved: de_all_results.RData\n"))

# Save a summary table across methods
summary_rows <- list()
for (cname in names(all_results)) {
  for (method in c("deseq2", "limma_voom")) {
    if (!is.null(all_results[[cname]][[method]])) {
      res_df <- all_results[[cname]][[method]]$results
      n_sig <- sum(!is.na(res_df$padj) & res_df$padj < opt$padj_threshold &
                     abs(res_df$log2FoldChange) > opt$lfc_threshold, na.rm = TRUE)
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        contrast = cname,
        method = method,
        n_tested = nrow(res_df),
        n_significant = n_sig,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  write.table(summary_df,
              file = file.path(opt$output_dir, "de_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("    Saved: de_summary.tsv\n")
  cat("\n    Summary:\n")
  print(summary_df)
}

cat("\n============================================================\n")
cat("Differential expression analysis complete!\n")
cat("============================================================\n")

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()),
           file.path(opt$output_dir, "session_info.txt"))
