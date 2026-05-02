#!/usr/bin/env Rscript
# ============================================================================
# PIPELINE VALIDATION CHECKS
# ============================================================================
# Performs 4 validation checks:
#   1. DE concordance: DESeq2 vs limma-voom agreement
#   2. Sex check: XIST / Y-gene expression vs reported sex
#   3. Expression outliers: PCA-based and correlation-based detection
#   4. Genomic inflation: lambda and QQ plots for each contrast
#
# Usage:
#   Rscript validation.R \
#     --deseq2_results_dir <path> \
#     --limma_results_dir <path> \
#     --count_matrix <path> \
#     --normalized_counts <path> \
#     --metadata <path> \
#     --output_dir <path> \
#     --outlier_sd_threshold <num> \
#     --min_sample_correlation <num> \
#     --threads <int>
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stats)
})

# ============================================================================
# PARSE ARGUMENTS
# ============================================================================

option_list <- list(
  make_option("--deseq2_results_dir", type = "character", default = NULL,
              help = "Directory containing DESeq2 result TSVs"),
  make_option("--limma_results_dir", type = "character", default = NULL,
              help = "Directory containing limma-voom result TSVs"),
  make_option("--count_matrix", type = "character", default = NULL,
              help = "Path to raw count matrix"),
  make_option("--normalized_counts", type = "character", default = NULL,
              help = "Path to normalized count matrix (VST or similar)"),
  make_option("--metadata", type = "character", default = NULL,
              help = "Path to sample metadata/samplesheet CSV"),
  make_option("--output_dir", type = "character", default = "validation_output",
              help = "Output directory"),
  make_option("--outlier_sd_threshold", type = "numeric", default = 3,
              help = "SD threshold for PCA outlier detection [default: 3]"),
  make_option("--min_sample_correlation", type = "numeric", default = 0.8,
              help = "Minimum mean pairwise correlation [default: 0.8]"),
  make_option("--threads", type = "integer", default = 4,
              help = "Number of threads [default: 4]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Create output subdirectories
plots_dir <- file.path(opt$output_dir, "plots")
reports_dir <- file.path(opt$output_dir, "reports")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

# Initialize report
report_lines <- character()
add_report <- function(...) {
  report_lines <<- c(report_lines, paste0(...))
}

add_report("# Pipeline Validation Report")
add_report("")
add_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
add_report("")

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Load normalized counts
norm_counts <- as.data.frame(fread(opt$normalized_counts), row.names = 1)
cat(sprintf("  Normalized counts: %d genes x %d samples\n", nrow(norm_counts), ncol(norm_counts)))

# Load raw counts
raw_counts <- as.data.frame(fread(opt$count_matrix), row.names = 1)
cat(sprintf("  Raw counts: %d genes x %d samples\n", nrow(raw_counts), ncol(raw_counts)))

# Load metadata
metadata <- fread(opt$metadata)
cat(sprintf("  Metadata: %d samples\n", nrow(metadata)))

# ============================================================================
# CHECK 1: DE CONCORDANCE (DESeq2 vs limma-voom)
# ============================================================================

add_report("## 1. Differential Expression Concordance")
add_report("")

run_de_concordance <- function(deseq2_dir, limma_dir) {
  if (!dir.exists(deseq2_dir) || !dir.exists(limma_dir)) {
    add_report("**SKIPPED**: DESeq2 or limma results directory not found.")
    add_report("")
    return(NULL)
  }

  deseq2_files <- list.files(deseq2_dir, pattern = "\\.tsv$", full.names = TRUE)
  limma_files <- list.files(limma_dir, pattern = "\\.tsv$", full.names = TRUE)

  if (length(deseq2_files) == 0 || length(limma_files) == 0) {
    add_report("**SKIPPED**: No result files found in DESeq2 or limma directories.")
    add_report("")
    return(NULL)
  }

  # Match contrasts by filename pattern
  deseq2_contrasts <- gsub("_deseq2\\.tsv$|\\.tsv$", "", basename(deseq2_files))
  limma_contrasts <- gsub("_limma\\.tsv$|\\.tsv$", "", basename(limma_files))

  common_contrasts <- intersect(deseq2_contrasts, limma_contrasts)

  if (length(common_contrasts) == 0) {
    # Try matching by position if names don't overlap
    n_match <- min(length(deseq2_files), length(limma_files))
    if (n_match > 0) {
      common_contrasts <- deseq2_contrasts[1:n_match]
      # Pair them by order
    } else {
      add_report("**WARNING**: Could not match any contrasts between DESeq2 and limma.")
      add_report("")
      return(NULL)
    }
  }

  concordance_results <- data.frame(
    contrast = character(),
    spearman_lfc = numeric(),
    jaccard_sig = numeric(),
    n_discordant = integer(),
    stringsAsFactors = FALSE
  )

  for (contrast in common_contrasts) {
    # Find matching files
    d_idx <- which(deseq2_contrasts == contrast)
    l_idx <- which(limma_contrasts == contrast)

    if (length(d_idx) == 0 || length(l_idx) == 0) next

    deseq2_res <- fread(deseq2_files[d_idx[1]])
    limma_res <- fread(limma_files[l_idx[1]])

    # Standardize column names
    # DESeq2 typically has: gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    # limma typically has: gene, logFC, AveExpr, t, P.Value, adj.P.Val, B
    deseq2_lfc_col <- grep("log2FoldChange|logFC|lfc", names(deseq2_res), value = TRUE, ignore.case = TRUE)[1]
    deseq2_padj_col <- grep("padj|adj.P.Val|FDR", names(deseq2_res), value = TRUE, ignore.case = TRUE)[1]
    limma_lfc_col <- grep("logFC|log2FoldChange|lfc", names(limma_res), value = TRUE, ignore.case = TRUE)[1]
    limma_padj_col <- grep("adj.P.Val|padj|FDR", names(limma_res), value = TRUE, ignore.case = TRUE)[1]

    # Get gene ID column
    gene_col_d <- grep("gene|Gene|gene_id|ensembl", names(deseq2_res), value = TRUE, ignore.case = TRUE)[1]
    gene_col_l <- grep("gene|Gene|gene_id|ensembl", names(limma_res), value = TRUE, ignore.case = TRUE)[1]

    if (is.na(gene_col_d)) gene_col_d <- names(deseq2_res)[1]
    if (is.na(gene_col_l)) gene_col_l <- names(limma_res)[1]

    if (any(is.na(c(deseq2_lfc_col, deseq2_padj_col, limma_lfc_col, limma_padj_col)))) {
      cat(sprintf("  WARNING: Cannot identify LFC/padj columns for contrast: %s\n", contrast))
      next
    }

    # Merge on gene
    merged <- merge(
      deseq2_res[, c(gene_col_d, deseq2_lfc_col, deseq2_padj_col), with = FALSE],
      limma_res[, c(gene_col_l, limma_lfc_col, limma_padj_col), with = FALSE],
      by.x = gene_col_d, by.y = gene_col_l
    )

    setnames(merged, c("gene", "deseq2_lfc", "deseq2_padj", "limma_lfc", "limma_padj"))

    # Remove NAs
    merged <- merged[complete.cases(merged)]

    if (nrow(merged) == 0) next

    # Spearman correlation of log2FC
    spearman_cor <- cor(merged$deseq2_lfc, merged$limma_lfc, method = "spearman")

    # Jaccard index of significant genes
    sig_deseq2 <- merged$gene[merged$deseq2_padj < 0.05]
    sig_limma <- merged$gene[merged$limma_padj < 0.05]
    jaccard <- length(intersect(sig_deseq2, sig_limma)) /
               length(union(sig_deseq2, sig_limma))

    # Discordant genes: significant in both but opposite direction
    both_sig <- intersect(sig_deseq2, sig_limma)
    merged_sig <- merged[merged$gene %in% both_sig, ]
    discordant <- sum(sign(merged_sig$deseq2_lfc) != sign(merged_sig$limma_lfc))

    concordance_results <- rbind(concordance_results, data.frame(
      contrast = contrast,
      spearman_lfc = round(spearman_cor, 4),
      jaccard_sig = round(jaccard, 4),
      n_discordant = discordant,
      stringsAsFactors = FALSE
    ))

    # Generate LFC comparison plot
    p <- ggplot(merged, aes(x = deseq2_lfc, y = limma_lfc)) +
      geom_point(alpha = 0.3, size = 0.5) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      labs(
        title = sprintf("DE Concordance: %s", contrast),
        subtitle = sprintf("Spearman rho = %.3f | Jaccard = %.3f | Discordant = %d",
                          spearman_cor, jaccard, discordant),
        x = "DESeq2 log2FC",
        y = "limma-voom logFC"
      ) +
      theme_minimal()

    ggsave(file.path(plots_dir, sprintf("concordance_%s.png", contrast)),
           p, width = 7, height = 6, dpi = 150)
  }

  # Report results
  if (nrow(concordance_results) > 0) {
    add_report("| Contrast | Spearman (LFC) | Jaccard (sig) | Discordant |")
    add_report("|----------|---------------|---------------|------------|")
    for (i in seq_len(nrow(concordance_results))) {
      row <- concordance_results[i, ]
      status <- ifelse(row$spearman_lfc < 0.7, " **LOW**", "")
      add_report(sprintf("| %s | %.3f%s | %.3f | %d |",
                        row$contrast, row$spearman_lfc, status, row$jaccard_sig, row$n_discordant))
    }
    add_report("")

    # Flag low concordance
    low_concordance <- concordance_results[concordance_results$spearman_lfc < 0.7, ]
    if (nrow(low_concordance) > 0) {
      add_report(sprintf("**WARNING**: %d contrast(s) have low DESeq2-limma concordance (rho < 0.7).",
                        nrow(low_concordance)))
      add_report("")
    }
  }

  return(concordance_results)
}

cat("Running DE concordance check...\n")
concordance <- run_de_concordance(opt$deseq2_results_dir, opt$limma_results_dir)

# ============================================================================
# CHECK 2: SEX CHECK
# ============================================================================

add_report("## 2. Sex Validation Check")
add_report("")

run_sex_check <- function(norm_counts, metadata) {
  # Sex-linked genes
  xist_gene <- "XIST"
  y_genes <- c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D")

  # Find sex column in metadata
  sex_col <- grep("^sex$|^gender$|^Sex$|^Gender$", names(metadata), value = TRUE)[1]

  if (is.na(sex_col)) {
    add_report("**SKIPPED**: No 'sex' or 'gender' column found in metadata.")
    add_report("")
    return(NULL)
  }

  # Get sample IDs from metadata
  sample_col <- names(metadata)[1]
  meta_samples <- metadata[[sample_col]]
  reported_sex <- metadata[[sex_col]]
  names(reported_sex) <- meta_samples

  # Find available genes in count matrix (try multiple naming conventions)
  available_genes <- rownames(norm_counts)

  find_gene <- function(gene_name) {
    # Try exact match
    idx <- grep(paste0("^", gene_name, "$"), available_genes, value = TRUE)
    if (length(idx) > 0) return(idx[1])
    # Try with ensembl suffix (e.g., XIST|ENSG...)
    idx <- grep(paste0("^", gene_name, "\\|"), available_genes, value = TRUE)
    if (length(idx) > 0) return(idx[1])
    idx <- grep(paste0("\\|", gene_name, "$"), available_genes, value = TRUE)
    if (length(idx) > 0) return(idx[1])
    # Try partial match
    idx <- grep(gene_name, available_genes, value = TRUE)
    if (length(idx) > 0) return(idx[1])
    return(NA)
  }

  xist_id <- find_gene(xist_gene)
  y_gene_ids <- sapply(y_genes, find_gene)
  y_gene_ids <- y_gene_ids[!is.na(y_gene_ids)]

  if (is.na(xist_id) && length(y_gene_ids) == 0) {
    add_report("**SKIPPED**: Neither XIST nor Y-chromosome genes found in count matrix.")
    add_report("")
    return(NULL)
  }

  # Common samples
  common_samples <- intersect(colnames(norm_counts), meta_samples)
  if (length(common_samples) == 0) {
    add_report("**SKIPPED**: No overlapping samples between counts and metadata.")
    add_report("")
    return(NULL)
  }

  # Calculate XIST expression
  sex_results <- data.frame(
    sample = common_samples,
    reported_sex = as.character(reported_sex[common_samples]),
    stringsAsFactors = FALSE
  )

  if (!is.na(xist_id)) {
    sex_results$xist_expr <- as.numeric(norm_counts[xist_id, common_samples])
  }

  # Calculate mean Y-gene expression
  if (length(y_gene_ids) > 0) {
    y_expr <- norm_counts[y_gene_ids, common_samples, drop = FALSE]
    sex_results$y_gene_mean <- colMeans(y_expr, na.rm = TRUE)
  }

  # Predict sex based on expression
  sex_results$predicted_sex <- NA
  if (!is.na(xist_id) && length(y_gene_ids) > 0) {
    xist_median <- median(sex_results$xist_expr, na.rm = TRUE)
    y_median <- median(sex_results$y_gene_mean, na.rm = TRUE)
    sex_results$predicted_sex <- ifelse(
      sex_results$xist_expr > xist_median & sex_results$y_gene_mean < y_median,
      "F", "M"
    )
  } else if (!is.na(xist_id)) {
    xist_median <- median(sex_results$xist_expr, na.rm = TRUE)
    sex_results$predicted_sex <- ifelse(sex_results$xist_expr > xist_median, "F", "M")
  } else {
    y_median <- median(sex_results$y_gene_mean, na.rm = TRUE)
    sex_results$predicted_sex <- ifelse(sex_results$y_gene_mean > y_median, "M", "F")
  }

  # Standardize reported sex for comparison
  sex_results$reported_std <- toupper(substr(sex_results$reported_sex, 1, 1))
  sex_results$reported_std[sex_results$reported_std == "MALE" | sex_results$reported_std == "M"] <- "M"
  sex_results$reported_std[sex_results$reported_std == "FEMALE" | sex_results$reported_std == "F"] <- "F"

  # Flag mismatches
  sex_results$mismatch <- (sex_results$predicted_sex != sex_results$reported_std) &
                          !is.na(sex_results$predicted_sex) &
                          sex_results$reported_std %in% c("M", "F")

  n_mismatches <- sum(sex_results$mismatch, na.rm = TRUE)
  n_checked <- sum(sex_results$reported_std %in% c("M", "F"))

  add_report(sprintf("- Samples checked: %d", n_checked))
  add_report(sprintf("- Sex mismatches detected: **%d**", n_mismatches))
  add_report("")

  if (n_mismatches > 0) {
    add_report("**MISMATCHED SAMPLES:**")
    add_report("")
    add_report("| Sample | Reported | Predicted | XIST | Y-gene mean |")
    add_report("|--------|----------|-----------|------|-------------|")
    mismatched <- sex_results[sex_results$mismatch == TRUE & !is.na(sex_results$mismatch), ]
    for (i in seq_len(nrow(mismatched))) {
      xist_val <- ifelse("xist_expr" %in% names(mismatched), sprintf("%.2f", mismatched$xist_expr[i]), "NA")
      y_val <- ifelse("y_gene_mean" %in% names(mismatched), sprintf("%.2f", mismatched$y_gene_mean[i]), "NA")
      add_report(sprintf("| %s | %s | %s | %s | %s |",
                        mismatched$sample[i], mismatched$reported_std[i],
                        mismatched$predicted_sex[i], xist_val, y_val))
    }
    add_report("")
    add_report("**ACTION REQUIRED**: Verify sample identity for mismatched samples.")
    add_report("")
  } else {
    add_report("All samples pass sex validation check.")
    add_report("")
  }

  # Generate sex check plot
  if (!is.na(xist_id) && length(y_gene_ids) > 0) {
    p <- ggplot(sex_results, aes(x = xist_expr, y = y_gene_mean,
                                  color = reported_std, shape = mismatch)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4),
                        labels = c("Match", "Mismatch")) +
      scale_color_manual(values = c("F" = "#E41A1C", "M" = "#377EB8"),
                        na.value = "grey50") +
      labs(
        title = "Sex Validation Check",
        x = "XIST Expression",
        y = "Mean Y-gene Expression",
        color = "Reported Sex",
        shape = "Status"
      ) +
      theme_minimal()

    ggsave(file.path(plots_dir, "sex_check.png"), p, width = 8, height = 6, dpi = 150)
  }

  return(sex_results)
}

cat("Running sex validation check...\n")
sex_check <- run_sex_check(norm_counts, metadata)

# ============================================================================
# CHECK 3: EXPRESSION OUTLIERS
# ============================================================================

add_report("## 3. Expression Outlier Detection")
add_report("")

run_outlier_detection <- function(norm_counts, outlier_sd, min_corr) {
  cat("Running PCA-based outlier detection...\n")

  # PCA on normalized counts (transpose: samples as rows)
  # Remove genes with zero variance
  gene_vars <- apply(norm_counts, 1, var, na.rm = TRUE)
  variable_genes <- names(gene_vars[gene_vars > 0])

  if (length(variable_genes) < 100) {
    add_report("**SKIPPED**: Fewer than 100 variable genes available for PCA.")
    add_report("")
    return(NULL)
  }

  # Use top 5000 most variable genes for efficiency
  top_genes <- names(sort(gene_vars[variable_genes], decreasing = TRUE))[1:min(5000, length(variable_genes))]
  pca_input <- t(norm_counts[top_genes, ])

  # Run PCA
  pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
  pca_scores <- as.data.frame(pca_result$x[, 1:min(10, ncol(pca_result$x))])
  pca_scores$sample <- rownames(pca_scores)

  # Calculate distance from centroid using first 5 PCs
  n_pcs <- min(5, ncol(pca_result$x))
  centroid <- colMeans(pca_scores[, 1:n_pcs])
  distances <- sqrt(rowSums(sweep(pca_scores[, 1:n_pcs], 2, centroid)^2))
  mean_dist <- mean(distances)
  sd_dist <- sd(distances)

  # Flag outliers
  pca_outliers <- names(distances[distances > mean_dist + outlier_sd * sd_dist])

  # Correlation-based outlier detection
  cat("Running correlation-based outlier detection...\n")
  sample_cor <- cor(norm_counts[top_genes, ], method = "spearman")
  mean_cor <- rowMeans(sample_cor) - 1/(ncol(sample_cor))  # subtract self-correlation contribution
  # More accurately: mean of off-diagonal
  diag(sample_cor) <- NA
  mean_cor_offdiag <- rowMeans(sample_cor, na.rm = TRUE)

  cor_outliers <- names(mean_cor_offdiag[mean_cor_offdiag < min_corr])

  # Combined outliers
  all_outliers <- union(pca_outliers, cor_outliers)

  add_report(sprintf("- Samples analyzed: %d", ncol(norm_counts)))
  add_report(sprintf("- PCA outliers (>%d SD from centroid): **%d**", outlier_sd, length(pca_outliers)))
  add_report(sprintf("- Correlation outliers (mean r < %.2f): **%d**", min_corr, length(cor_outliers)))
  add_report(sprintf("- Total unique outliers: **%d**", length(all_outliers)))
  add_report("")


  if (length(all_outliers) > 0) {
    add_report("**FLAGGED SAMPLES:**")
    add_report("")
    add_report("| Sample | PCA Distance (SD) | Mean Correlation | Flags |")
    add_report("|--------|-------------------|-----------------|-------|")
    for (s in all_outliers) {
      dist_sd <- (distances[s] - mean_dist) / sd_dist
      mcor <- mean_cor_offdiag[s]
      flags <- c()
      if (s %in% pca_outliers) flags <- c(flags, "PCA")
      if (s %in% cor_outliers) flags <- c(flags, "Correlation")
      add_report(sprintf("| %s | %.2f | %.3f | %s |", s, dist_sd, mcor, paste(flags, collapse = ", ")))
    }
    add_report("")
  } else {
    add_report("No expression outliers detected.")
    add_report("")
  }

  # Generate PCA outlier plot
  pca_scores$outlier <- pca_scores$sample %in% all_outliers
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100

  p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = outlier)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red"),
                      labels = c("Normal", "Outlier")) +
    labs(
      title = "Expression PCA - Outlier Detection",
      subtitle = sprintf("%d outlier(s) flagged", length(all_outliers)),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      color = "Status"
    ) +
    theme_minimal()

  # Label outliers
  if (length(all_outliers) > 0 && length(all_outliers) <= 20) {
    outlier_df <- pca_scores[pca_scores$outlier, ]
    p <- p + geom_text(data = outlier_df, aes(label = sample),
                       vjust = -0.5, size = 2.5, color = "red")
  }

  ggsave(file.path(plots_dir, "pca_outliers.png"), p, width = 9, height = 7, dpi = 150)

  # Correlation heatmap data for report
  return(list(pca_outliers = pca_outliers, cor_outliers = cor_outliers,
              all_outliers = all_outliers))
}

outliers <- run_outlier_detection(norm_counts, opt$outlier_sd_threshold, opt$min_sample_correlation)

# ============================================================================
# CHECK 4: GENOMIC INFLATION
# ============================================================================

add_report("## 4. Genomic Inflation (Lambda)")
add_report("")

run_genomic_inflation <- function(deseq2_dir) {
  if (!dir.exists(deseq2_dir)) {
    add_report("**SKIPPED**: DESeq2 results directory not found.")
    add_report("")
    return(NULL)
  }

  result_files <- list.files(deseq2_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(result_files) == 0) {
    add_report("**SKIPPED**: No DESeq2 result files found.")
    add_report("")
    return(NULL)
  }

  lambda_results <- data.frame(
    contrast = character(),
    lambda = numeric(),
    n_genes = integer(),
    status = character(),
    stringsAsFactors = FALSE
  )

  for (f in result_files) {
    contrast_name <- gsub("_deseq2\\.tsv$|\\.tsv$", "", basename(f))
    res <- fread(f)

    # Find pvalue column
    pval_col <- grep("^pvalue$|^P.Value$|^p.value$", names(res), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(pval_col)) {
      pval_col <- grep("pval", names(res), value = TRUE, ignore.case = TRUE)[1]
    }

    if (is.na(pval_col)) {
      cat(sprintf("  WARNING: No p-value column found in %s\n", basename(f)))
      next
    }

    pvals <- res[[pval_col]]
    pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals < 1]

    if (length(pvals) < 100) next

    # Calculate lambda (genomic inflation factor)
    # lambda = median(chi-squared statistics) / expected median of chi-sq(1)
    chisq_obs <- qchisq(1 - pvals, df = 1)
    lambda <- median(chisq_obs) / qchisq(0.5, df = 1)

    status <- "PASS"
    if (lambda > 1.5) status <- "FAIL (inflated)"
    if (lambda < 0.5) status <- "FAIL (deflated)"
    if (lambda > 1.2 && lambda <= 1.5) status <- "WARNING"

    lambda_results <- rbind(lambda_results, data.frame(
      contrast = contrast_name,
      lambda = round(lambda, 4),
      n_genes = length(pvals),
      status = status,
      stringsAsFactors = FALSE
    ))

    # Generate QQ plot
    n <- length(pvals)
    expected <- -log10(ppoints(n))
    observed <- sort(-log10(pvals))

    qq_df <- data.frame(
      expected = expected,
      observed = observed
    )

    p <- ggplot(qq_df, aes(x = expected, y = observed)) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      labs(
        title = sprintf("QQ Plot: %s", contrast_name),
        subtitle = sprintf("lambda = %.3f (%s)", lambda, status),
        x = expression(-log[10](expected~p)),
        y = expression(-log[10](observed~p))
      ) +
      theme_minimal() +
      coord_equal(ratio = max(expected) / max(observed))

    ggsave(file.path(plots_dir, sprintf("qq_%s.png", contrast_name)),
           p, width = 6, height = 6, dpi = 150)
  }

  # Report
  if (nrow(lambda_results) > 0) {
    add_report("| Contrast | Lambda | N genes | Status |")
    add_report("|----------|--------|---------|--------|")
    for (i in seq_len(nrow(lambda_results))) {
      row <- lambda_results[i, ]
      add_report(sprintf("| %s | %.3f | %d | %s |", row$contrast, row$lambda, row$n_genes, row$status))
    }
    add_report("")

    inflated <- lambda_results[grepl("FAIL|WARNING", lambda_results$status), ]
    if (nrow(inflated) > 0) {
      add_report(sprintf("**WARNING**: %d contrast(s) show genomic inflation/deflation.", nrow(inflated)))
      add_report("Consider reviewing covariates or batch correction.")
      add_report("")
    } else {
      add_report("All contrasts show acceptable genomic inflation (lambda between 0.5 and 1.5).")
      add_report("")
    }
  }

  return(lambda_results)
}

cat("Running genomic inflation analysis...\n")
inflation <- run_genomic_inflation(opt$deseq2_results_dir)

# ============================================================================
# COMPILE FINAL REPORT
# ============================================================================

add_report("---")
add_report("")
add_report("## Summary")
add_report("")

# Count issues
n_issues <- 0
if (!is.null(concordance) && any(concordance$spearman_lfc < 0.7)) {
  n_issues <- n_issues + sum(concordance$spearman_lfc < 0.7)
}
if (!is.null(sex_check) && any(sex_check$mismatch, na.rm = TRUE)) {
  n_issues <- n_issues + sum(sex_check$mismatch, na.rm = TRUE)
}
if (!is.null(outliers)) {
  n_issues <- n_issues + length(outliers$all_outliers)
}
if (!is.null(inflation) && any(grepl("FAIL", inflation$status))) {
  n_issues <- n_issues + sum(grepl("FAIL", inflation$status))
}

if (n_issues == 0) {
  add_report("**PASS**: All validation checks passed with no critical issues.")
} else {
  add_report(sprintf("**ATTENTION**: %d issue(s) flagged across all checks. Review above sections for details.", n_issues))
}
add_report("")
add_report(sprintf("Plots saved to: `%s`", plots_dir))
add_report("")

# Write report
report_file <- file.path(reports_dir, "validation_report.md")
writeLines(report_lines, report_file)
cat(sprintf("\nValidation report written to: %s\n", report_file))

# Also save results as RData for programmatic access
save(concordance, sex_check, outliers, inflation,
     file = file.path(opt$output_dir, "validation_results.RData"))

cat("\nValidation complete.\n")
