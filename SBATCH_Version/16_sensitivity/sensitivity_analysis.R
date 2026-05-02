#!/usr/bin/env Rscript
# ============================================================================
# 16 - Sensitivity Analysis
# ============================================================================
# Evaluates robustness of differential expression results via:
# 1. Timepoint-stratified DE (diagnostic-only, relapse-only subsets)
# 2. Ancestry threshold sweep (vary proportion cutoffs, measure DE stability)
# 3. Covariate leave-one-out (remove each covariate, re-run DE, Jaccard overlap)
# 4. Bootstrap DE stability (resample, measure gene list stability)
#
# Required packages: DESeq2, optparse, tidyverse, UpSetR, parallel
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(tidyverse)
    library(UpSetR)
    library(parallel)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
    make_option("--count_matrix", type = "character",
                help = "Path to raw count matrix (genes x samples TSV)"),
    make_option("--metadata", type = "character",
                help = "Path to sample metadata CSV"),
    make_option("--ancestry_proportions", type = "character",
                help = "Path to ancestry proportions TSV"),
    make_option("--de_results_dir", type = "character",
                help = "Directory containing original DE results"),
    make_option("--output_dir", type = "character",
                help = "Output directory"),
    make_option("--bootstrap_iterations", type = "integer", default = 1000,
                help = "Number of bootstrap iterations [default: %default]"),
    make_option("--covariates", type = "character", default = "batch,sex,age",
                help = "Comma-separated covariates used in DE [default: %default]"),
    make_option("--threads", type = "integer", default = 8,
                help = "Number of threads [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$count_matrix) || is.null(opt$metadata) ||
    is.null(opt$ancestry_proportions) || is.null(opt$de_results_dir) ||
    is.null(opt$output_dir)) {
    stop("Required: --count_matrix, --metadata, --ancestry_proportions, --de_results_dir, --output_dir")
}

covariates <- strsplit(opt$covariates, ",")[[1]]
n_boot <- opt$bootstrap_iterations

cat("=== Sensitivity Analysis ===\n")
cat("Count matrix:", opt$count_matrix, "\n")
cat("Metadata:", opt$metadata, "\n")
cat("Ancestry proportions:", opt$ancestry_proportions, "\n")
cat("DE results dir:", opt$de_results_dir, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Bootstrap iterations:", n_boot, "\n")
cat("Covariates:", paste(covariates, collapse = ", "), "\n")
cat("Threads:", opt$threads, "\n\n")

# ============================================================================
# Load data
# ============================================================================
cat("Loading data...\n")

counts <- read.table(opt$count_matrix, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE)
metadata <- read.csv(opt$metadata, header = TRUE, stringsAsFactors = FALSE)
ancestry <- read.table(opt$ancestry_proportions, header = TRUE, sep = "\t",
                       row.names = 1, check.names = FALSE)

cat(sprintf("  Count matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  Metadata: %d samples\n", nrow(metadata)))
cat(sprintf("  Ancestry proportions: %d samples x %d components\n", nrow(ancestry), ncol(ancestry)))

# Load original DE results for comparison
original_de_files <- list.files(opt$de_results_dir, pattern = "*.tsv",
                                full.names = TRUE, recursive = TRUE)
cat(sprintf("  Found %d original DE result files\n", length(original_de_files)))

# Load first DE result as reference
if (length(original_de_files) > 0) {
    original_de <- read.table(original_de_files[1], header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE)
    original_sig_genes <- original_de %>%
        filter(padj < 0.05, abs(log2FoldChange) > 0.585) %>%
        pull(1)  # Gene names in first column
    cat(sprintf("  Reference DE: %d significant genes\n", length(original_sig_genes)))
} else {
    original_sig_genes <- character(0)
    cat("  WARNING: No original DE results found for comparison\n")
}

# Align metadata and counts
sample_ids <- colnames(counts)
metadata <- metadata %>% filter(metadata[[1]] %in% sample_ids)
rownames(metadata) <- metadata[[1]]
metadata <- metadata[sample_ids[sample_ids %in% rownames(metadata)], ]
counts <- counts[, rownames(metadata)]

# Merge ancestry into metadata
if (any(rownames(ancestry) %in% rownames(metadata))) {
    ancestry_aligned <- ancestry[rownames(metadata), , drop = FALSE]
    metadata <- cbind(metadata, ancestry_aligned)
}

# ============================================================================
# Helper: Run DESeq2 on a subset
# ============================================================================
run_deseq2_subset <- function(count_mat, meta, design_formula, padj_cut = 0.05, lfc_cut = 0.585) {
    tryCatch({
        dds <- DESeqDataSetFromMatrix(
            countData = count_mat,
            colData = meta,
            design = design_formula
        )
        # Filter low-count genes
        keep <- rowSums(counts(dds) >= 10) >= 3
        dds <- dds[keep, ]
        dds <- DESeq(dds, parallel = FALSE, quiet = TRUE)
        res <- results(dds, alpha = padj_cut)
        sig_genes <- rownames(res)[which(res$padj < padj_cut & abs(res$log2FoldChange) > lfc_cut)]
        return(sig_genes)
    }, error = function(e) {
        cat(sprintf("    DESeq2 subset failed: %s\n", conditionMessage(e)))
        return(character(0))
    })
}

# Helper: Compute Jaccard index
jaccard <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union_size <- length(union(set1, set2))
    if (union_size == 0) return(0)
    return(intersection / union_size)
}

# ============================================================================
# 1. Timepoint-stratified DE
# ============================================================================
cat("\n--- 1. Timepoint-Stratified DE ---\n")
timepoint_results <- list()

# Identify timepoint column
timepoint_col <- NULL
for (col in c("timepoint", "Timepoint", "time_point", "diagnosis_type")) {
    if (col %in% colnames(metadata)) {
        timepoint_col <- col
        break
    }
}

if (!is.null(timepoint_col)) {
    timepoints <- unique(metadata[[timepoint_col]])
    cat(sprintf("  Timepoint column: %s\n", timepoint_col))
    cat(sprintf("  Timepoints found: %s\n", paste(timepoints, collapse = ", ")))

    # Identify condition column (first column after sample_id that looks like group/condition)
    condition_col <- NULL
    for (col in c("condition", "group", "diagnosis", "subtype")) {
        if (col %in% colnames(metadata)) {
            condition_col <- col
            break
        }
    }

    if (!is.null(condition_col) && length(unique(metadata[[condition_col]])) >= 2) {
        for (tp in timepoints) {
            cat(sprintf("  Running DE for timepoint: %s\n", tp))
            tp_mask <- metadata[[timepoint_col]] == tp
            tp_meta <- metadata[tp_mask, , drop = FALSE]
            tp_counts <- counts[, rownames(tp_meta), drop = FALSE]

            if (nrow(tp_meta) >= 6 && length(unique(tp_meta[[condition_col]])) >= 2) {
                design_formula <- as.formula(paste("~", condition_col))
                sig_genes <- run_deseq2_subset(tp_counts, tp_meta, design_formula)
                timepoint_results[[tp]] <- sig_genes
                cat(sprintf("    %s: %d significant genes\n", tp, length(sig_genes)))

                # Save results
                write.table(data.frame(gene = sig_genes),
                            file.path(opt$output_dir, "timepoint_stratified",
                                      paste0("sig_genes_", tp, ".tsv")),
                            sep = "\t", quote = FALSE, row.names = FALSE)
            } else {
                cat(sprintf("    %s: insufficient samples or groups, skipping\n", tp))
            }
        }
    } else {
        cat("  WARNING: No suitable condition column found for stratified DE\n")
    }
} else {
    cat("  WARNING: No timepoint column found in metadata\n")
}

# Jaccard overlap between timepoints and original
if (length(timepoint_results) > 0 && length(original_sig_genes) > 0) {
    tp_jaccard <- data.frame(
        timepoint = names(timepoint_results),
        n_sig_genes = sapply(timepoint_results, length),
        jaccard_vs_original = sapply(timepoint_results, function(x) jaccard(x, original_sig_genes))
    )
    write.table(tp_jaccard,
                file.path(opt$output_dir, "timepoint_stratified", "timepoint_jaccard_overlap.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

# ============================================================================
# 2. Ancestry Threshold Sweep
# ============================================================================
cat("\n--- 2. Ancestry Threshold Sweep ---\n")

ancestry_cols <- colnames(ancestry)
sweep_results <- list()

# Define thresholds to sweep (proportion cutoffs for majority ancestry assignment)
thresholds <- seq(0.5, 0.9, by = 0.05)
cat(sprintf("  Ancestry components: %s\n", paste(ancestry_cols, collapse = ", ")))
cat(sprintf("  Thresholds: %s\n", paste(thresholds, collapse = ", ")))

# For each ancestry component, sweep thresholds
for (anc_col in ancestry_cols[1:min(3, length(ancestry_cols))]) {
    cat(sprintf("  Sweeping thresholds for: %s\n", anc_col))

    for (thresh in thresholds) {
        # Assign samples to high/low ancestry groups based on threshold
        if (anc_col %in% colnames(metadata)) {
            high_anc <- metadata[[anc_col]] >= thresh
            low_anc <- metadata[[anc_col]] < (1 - thresh)

            # Only run if we have enough samples in both groups
            if (sum(high_anc) >= 3 && sum(low_anc) >= 3) {
                subset_mask <- high_anc | low_anc
                sub_meta <- metadata[subset_mask, , drop = FALSE]
                sub_meta$ancestry_group <- ifelse(
                    sub_meta[[anc_col]] >= thresh, "high", "low"
                )
                sub_counts <- counts[, rownames(sub_meta), drop = FALSE]

                design_formula <- as.formula("~ ancestry_group")
                sig_genes <- run_deseq2_subset(sub_counts, sub_meta, design_formula)

                sweep_results[[paste(anc_col, thresh, sep = "_")]] <- list(
                    ancestry = anc_col,
                    threshold = thresh,
                    n_high = sum(high_anc),
                    n_low = sum(low_anc),
                    n_sig = length(sig_genes),
                    jaccard = jaccard(sig_genes, original_sig_genes),
                    genes = sig_genes
                )
            }
        }
    }
}

# Save sweep summary
if (length(sweep_results) > 0) {
    sweep_summary <- do.call(rbind, lapply(sweep_results, function(x) {
        data.frame(ancestry = x$ancestry, threshold = x$threshold,
                   n_high = x$n_high, n_low = x$n_low,
                   n_significant_genes = x$n_sig,
                   jaccard_vs_original = x$jaccard,
                   stringsAsFactors = FALSE)
    }))
    write.table(sweep_summary,
                file.path(opt$output_dir, "ancestry_sweep", "ancestry_threshold_sweep.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Sweep results saved (%d configurations tested)\n", nrow(sweep_summary)))
}

# ============================================================================
# 3. Covariate Leave-One-Out
# ============================================================================
cat("\n--- 3. Covariate Leave-One-Out ---\n")

loo_results <- list()

# Identify condition column
condition_col_loo <- NULL
for (col in c("condition", "group", "diagnosis", "subtype")) {
    if (col %in% colnames(metadata)) {
        condition_col_loo <- col
        break
    }
}

if (!is.null(condition_col_loo)) {
    # Full model with all covariates
    available_covs <- covariates[covariates %in% colnames(metadata)]
    cat(sprintf("  Condition: %s\n", condition_col_loo))
    cat(sprintf("  Available covariates: %s\n", paste(available_covs, collapse = ", ")))

    # Run full model
    if (length(available_covs) > 0) {
        full_formula <- as.formula(paste("~", paste(c(available_covs, condition_col_loo), collapse = " + ")))
        cat("  Running full model...\n")
        full_sig <- run_deseq2_subset(counts[, rownames(metadata)], metadata, full_formula)
        loo_results[["full_model"]] <- full_sig
        cat(sprintf("    Full model: %d significant genes\n", length(full_sig)))

        # Leave-one-out for each covariate
        for (cov in available_covs) {
            cat(sprintf("  Removing covariate: %s\n", cov))
            reduced_covs <- available_covs[available_covs != cov]

            if (length(reduced_covs) > 0) {
                reduced_formula <- as.formula(paste("~", paste(c(reduced_covs, condition_col_loo), collapse = " + ")))
            } else {
                reduced_formula <- as.formula(paste("~", condition_col_loo))
            }

            loo_sig <- run_deseq2_subset(counts[, rownames(metadata)], metadata, reduced_formula)
            loo_results[[paste0("without_", cov)]] <- loo_sig
            cat(sprintf("    Without %s: %d sig genes (Jaccard vs full: %.3f)\n",
                        cov, length(loo_sig), jaccard(loo_sig, full_sig)))
        }

        # Save LOO summary
        loo_summary <- data.frame(
            model = names(loo_results),
            n_significant_genes = sapply(loo_results, length),
            jaccard_vs_full = sapply(loo_results, function(x) jaccard(x, full_sig)),
            jaccard_vs_original = sapply(loo_results, function(x) jaccard(x, original_sig_genes)),
            stringsAsFactors = FALSE
        )
        write.table(loo_summary,
                    file.path(opt$output_dir, "leave_one_out", "loo_summary.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
        cat("  WARNING: No covariates available in metadata, skipping LOO\n")
    }
} else {
    cat("  WARNING: No condition column found, skipping LOO\n")
}

# ============================================================================
# 4. Bootstrap DE Stability
# ============================================================================
cat("\n--- 4. Bootstrap DE Stability ---\n")
cat(sprintf("  Running %d bootstrap iterations...\n", n_boot))

# Use parallel processing
n_cores <- min(opt$threads, detectCores())
cat(sprintf("  Using %d cores\n", n_cores))

if (!is.null(condition_col_loo)) {
    # Bootstrap: resample with replacement, run DE, track gene selection frequency
    set.seed(42)

    bootstrap_one <- function(iter) {
        # Resample samples with replacement within each group
        groups <- unique(metadata[[condition_col_loo]])
        boot_indices <- c()
        for (grp in groups) {
            grp_idx <- which(metadata[[condition_col_loo]] == grp)
            boot_indices <- c(boot_indices, sample(grp_idx, length(grp_idx), replace = TRUE))
        }

        # Create bootstrap dataset (handle duplicate indices)
        boot_meta <- metadata[boot_indices, , drop = FALSE]
        boot_counts <- counts[, boot_indices, drop = FALSE]

        # Make unique sample names for duplicates
        rownames(boot_meta) <- paste0("s", seq_len(nrow(boot_meta)))
        colnames(boot_counts) <- rownames(boot_meta)

        design_formula <- as.formula(paste("~", condition_col_loo))
        sig_genes <- run_deseq2_subset(boot_counts, boot_meta, design_formula)
        return(sig_genes)
    }

    # Run bootstrap iterations (with progress reporting)
    boot_results <- mclapply(seq_len(n_boot), function(i) {
        if (i %% 100 == 0) cat(sprintf("    Iteration %d / %d\n", i, n_boot))
        bootstrap_one(i)
    }, mc.cores = n_cores)

    # Compute gene selection frequency
    all_boot_genes <- unlist(boot_results)
    gene_freq <- sort(table(all_boot_genes) / n_boot, decreasing = TRUE)

    boot_stability <- data.frame(
        gene = names(gene_freq),
        selection_frequency = as.numeric(gene_freq),
        stringsAsFactors = FALSE
    )
    boot_stability$in_original <- boot_stability$gene %in% original_sig_genes

    write.table(boot_stability,
                file.path(opt$output_dir, "bootstrap", "gene_selection_frequency.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Compute pairwise Jaccard between bootstrap iterations (subsample for speed)
    n_compare <- min(100, n_boot)
    boot_jaccards <- numeric(n_compare * (n_compare - 1) / 2)
    idx <- 1
    for (i in 1:(n_compare - 1)) {
        for (j in (i + 1):n_compare) {
            boot_jaccards[idx] <- jaccard(boot_results[[i]], boot_results[[j]])
            idx <- idx + 1
        }
    }

    boot_metrics <- data.frame(
        metric = c("mean_jaccard", "median_jaccard", "sd_jaccard",
                   "mean_n_genes", "median_n_genes",
                   "genes_selected_gt50pct", "genes_selected_gt80pct",
                   "genes_selected_gt95pct"),
        value = c(mean(boot_jaccards), median(boot_jaccards), sd(boot_jaccards),
                  mean(sapply(boot_results, length)),
                  median(sapply(boot_results, length)),
                  sum(gene_freq > 0.5), sum(gene_freq > 0.8), sum(gene_freq > 0.95))
    )
    write.table(boot_metrics,
                file.path(opt$output_dir, "bootstrap", "bootstrap_stability_metrics.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    cat(sprintf("  Mean pairwise Jaccard: %.3f\n", mean(boot_jaccards)))
    cat(sprintf("  Genes selected in >50%% bootstraps: %d\n", sum(gene_freq > 0.5)))
    cat(sprintf("  Genes selected in >80%% bootstraps: %d\n", sum(gene_freq > 0.8)))
} else {
    cat("  WARNING: No condition column found, skipping bootstrap\n")
    boot_results <- list()
    loo_results <- list()
}

# ============================================================================
# Generate Plots
# ============================================================================
cat("\n--- Generating Plots ---\n")
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --- UpSet plot of gene overlap across analyses ---
if (length(timepoint_results) > 0 || length(loo_results) > 0) {
    cat("  Creating UpSet plot of gene overlaps...\n")
    tryCatch({
        # Combine all gene sets for UpSet
        all_gene_sets <- c(
            list(original = original_sig_genes),
            timepoint_results,
            loo_results
        )
        # Filter to non-empty sets
        all_gene_sets <- all_gene_sets[sapply(all_gene_sets, length) > 0]

        if (length(all_gene_sets) >= 2) {
            # Create binary matrix for UpSet
            all_genes <- unique(unlist(all_gene_sets))
            upset_mat <- data.frame(matrix(0, nrow = length(all_genes), ncol = length(all_gene_sets)))
            rownames(upset_mat) <- all_genes
            colnames(upset_mat) <- names(all_gene_sets)

            for (s in names(all_gene_sets)) {
                upset_mat[all_gene_sets[[s]], s] <- 1
            }

            pdf(file.path(plot_dir, "upset_gene_overlap.pdf"), width = 12, height = 8)
            print(upset(upset_mat, sets = colnames(upset_mat),
                        order.by = "freq", nsets = length(all_gene_sets),
                        mainbar.y.label = "Intersection Size",
                        sets.x.label = "Set Size"))
            dev.off()
            cat("    UpSet plot saved.\n")
        }
    }, error = function(e) {
        cat(sprintf("    WARNING: UpSet plot failed: %s\n", conditionMessage(e)))
    })
}

# --- Ancestry threshold sweep plot ---
if (length(sweep_results) > 0) {
    cat("  Creating ancestry sweep plot...\n")
    tryCatch({
        p_sweep <- ggplot(sweep_summary, aes(x = threshold, y = n_significant_genes,
                                              color = ancestry)) +
            geom_line(linewidth = 1) +
            geom_point(size = 2) +
            theme_minimal() +
            labs(
                title = "Ancestry Threshold Sweep: DE Gene Count Stability",
                x = "Ancestry Proportion Threshold",
                y = "Number of Significant DE Genes",
                color = "Ancestry Component"
            )
        ggsave(file.path(plot_dir, "ancestry_threshold_sweep.pdf"),
               p_sweep, width = 10, height = 6)

        # Jaccard stability plot
        p_jaccard <- ggplot(sweep_summary, aes(x = threshold, y = jaccard_vs_original,
                                                color = ancestry)) +
            geom_line(linewidth = 1) +
            geom_point(size = 2) +
            theme_minimal() +
            ylim(0, 1) +
            labs(
                title = "Ancestry Threshold Sweep: Jaccard Similarity vs Original DE",
                x = "Ancestry Proportion Threshold",
                y = "Jaccard Similarity",
                color = "Ancestry Component"
            )
        ggsave(file.path(plot_dir, "ancestry_sweep_jaccard.pdf"),
               p_jaccard, width = 10, height = 6)
        cat("    Ancestry sweep plots saved.\n")
    }, error = function(e) {
        cat(sprintf("    WARNING: Sweep plot failed: %s\n", conditionMessage(e)))
    })
}

# --- Bootstrap stability plot ---
if (exists("boot_stability") && nrow(boot_stability) > 0) {
    cat("  Creating bootstrap stability plot...\n")
    tryCatch({
        top_genes <- boot_stability %>% slice_head(n = 50)
        top_genes$gene <- factor(top_genes$gene, levels = rev(top_genes$gene))

        p_boot <- ggplot(top_genes, aes(x = selection_frequency, y = gene,
                                         fill = in_original)) +
            geom_bar(stat = "identity") +
            theme_minimal() +
            scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "coral")) +
            labs(
                title = "Top 50 Genes by Bootstrap Selection Frequency",
                x = "Selection Frequency",
                y = "Gene",
                fill = "In Original DE"
            )
        ggsave(file.path(plot_dir, "bootstrap_gene_stability.pdf"),
               p_boot, width = 10, height = 12)
        cat("    Bootstrap stability plot saved.\n")
    }, error = function(e) {
        cat(sprintf("    WARNING: Bootstrap plot failed: %s\n", conditionMessage(e)))
    })
}

# --- LOO Jaccard barplot ---
if (exists("loo_summary") && nrow(loo_summary) > 0) {
    cat("  Creating LOO Jaccard plot...\n")
    tryCatch({
        loo_plot_data <- loo_summary %>% filter(model != "full_model")
        loo_plot_data$model <- gsub("without_", "", loo_plot_data$model)

        p_loo <- ggplot(loo_plot_data, aes(x = reorder(model, -jaccard_vs_full),
                                            y = jaccard_vs_full)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ylim(0, 1) +
            labs(
                title = "Covariate Leave-One-Out: Jaccard Similarity vs Full Model",
                x = "Removed Covariate",
                y = "Jaccard Similarity"
            )
        ggsave(file.path(plot_dir, "loo_jaccard_barplot.pdf"),
               p_loo, width = 8, height = 6)
        cat("    LOO Jaccard plot saved.\n")
    }, error = function(e) {
        cat(sprintf("    WARNING: LOO plot failed: %s\n", conditionMessage(e)))
    })
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Sensitivity Analysis Complete ===\n")
cat(sprintf("Output directory: %s\n", opt$output_dir))
cat("Key outputs:\n")
cat(sprintf("  - Timepoint stratified: %s/timepoint_stratified/\n", opt$output_dir))
cat(sprintf("  - Ancestry sweep: %s/ancestry_sweep/\n", opt$output_dir))
cat(sprintf("  - Leave-one-out: %s/leave_one_out/\n", opt$output_dir))
cat(sprintf("  - Bootstrap: %s/bootstrap/\n", opt$output_dir))
cat(sprintf("  - Plots: %s/plots/\n", opt$output_dir))
