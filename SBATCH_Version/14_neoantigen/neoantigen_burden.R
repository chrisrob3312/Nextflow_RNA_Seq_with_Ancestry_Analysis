#!/usr/bin/env Rscript
# ============================================================================
# 14 - Neoantigen Burden Analysis
# ============================================================================
# Calculates per-sample neoantigen burden stratified by source (SNV, fusion,
# splicing, HERV) and compares across clinical and ancestry groups:
#   - Relapse vs non-relapse (Wilcoxon rank-sum)
#   - MRD-positive vs MRD-negative (Wilcoxon rank-sum)
#   - Ancestry groups (Kruskal-Wallis)
# Generates boxplots, stacked bar charts, and an annotated burden table.
#
# Required packages: ggplot2, tidyverse, optparse
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(tidyverse)
    library(ggplot2)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
    make_option("--merged_neoantigens", type = "character",
                help = "Path to merged neoantigen summary TSV (sample, pvacseq_count, neofuse_count, snaf_count, total)"),
    make_option("--metadata", type = "character",
                help = "Path to sample metadata CSV"),
    make_option("--ancestry_proportions", type = "character", default = "none",
                help = "Path to ancestry proportions TSV [default: none]"),
    make_option("--output_dir", type = "character",
                help = "Output directory"),
    make_option("--herv_dir", type = "character", default = "none",
                help = "Path to HERV quantification directory [default: none]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$merged_neoantigens) || is.null(opt$metadata) ||
    is.null(opt$output_dir)) {
    stop("Required: --merged_neoantigens, --metadata, --output_dir")
}

cat("=== Neoantigen Burden Analysis ===\n")
cat("Merged neoantigens:", opt$merged_neoantigens, "\n")
cat("Metadata:", opt$metadata, "\n")
cat("Ancestry proportions:", opt$ancestry_proportions, "\n")
cat("Output directory:", opt$output_dir, "\n\n")

# ============================================================================
# Create output directories
# ============================================================================
comp_dir <- file.path(opt$output_dir, "comparisons")
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Load data
# ============================================================================
cat("Loading data...\n")

# Load neoantigen summary
burden_df <- read.table(opt$merged_neoantigens, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
cat(sprintf("  Neoantigen summary: %d samples\n", nrow(burden_df)))

# Standardize column names
sample_col <- colnames(burden_df)[1]
colnames(burden_df)[1] <- "sample"

# Rename source columns to standardized names if needed
source_map <- c(
    pvacseq_count = "SNV",
    neofuse_count = "Fusion",
    snaf_count = "Splicing"
)
for (old_name in names(source_map)) {
    if (old_name %in% colnames(burden_df)) {
        colnames(burden_df)[colnames(burden_df) == old_name] <- source_map[old_name]
    }
}

# Add HERV burden if available
if (opt$herv_dir != "none" && dir.exists(opt$herv_dir)) {
    cat("  Loading HERV counts...\n")
    herv_counts <- data.frame(sample = character(0), HERV = integer(0))
    for (sid in burden_df$sample) {
        herv_file <- file.path(opt$herv_dir, sid, paste0(sid, "_herv_counts.tsv"))
        if (file.exists(herv_file)) {
            n_herv <- nrow(read.table(herv_file, header = FALSE, sep = "\t",
                                      stringsAsFactors = FALSE, fill = TRUE))
            herv_counts <- rbind(herv_counts, data.frame(sample = sid, HERV = n_herv))
        } else {
            herv_counts <- rbind(herv_counts, data.frame(sample = sid, HERV = 0L))
        }
    }
    burden_df <- burden_df %>% left_join(herv_counts, by = "sample")
} else {
    burden_df$HERV <- 0L
}

# Recalculate total if source columns exist
source_cols <- intersect(c("SNV", "Fusion", "Splicing", "HERV"), colnames(burden_df))
if (length(source_cols) > 0) {
    burden_df$total_burden <- rowSums(burden_df[, source_cols, drop = FALSE],
                                      na.rm = TRUE)
} else if ("total_neoantigens" %in% colnames(burden_df)) {
    burden_df$total_burden <- burden_df$total_neoantigens
} else {
    stop("Cannot determine neoantigen burden: no source or total column found")
}

cat(sprintf("  Burden sources: %s\n", paste(source_cols, collapse = ", ")))

# Load metadata
metadata <- read.csv(opt$metadata, header = TRUE, stringsAsFactors = FALSE)
meta_sample_col <- colnames(metadata)[1]
colnames(metadata)[1] <- "sample"

# Merge burden with metadata
burden_df <- burden_df %>%
    left_join(metadata, by = "sample")

cat(sprintf("  Merged data: %d samples\n", nrow(burden_df)))

# Load ancestry proportions
ancestry <- NULL
majority_ancestry <- NULL
if (opt$ancestry_proportions != "none" && file.exists(opt$ancestry_proportions)) {
    ancestry <- read.table(opt$ancestry_proportions, header = TRUE, sep = "\t",
                           row.names = 1, check.names = FALSE)
    common_samples <- intersect(rownames(ancestry), burden_df$sample)
    if (length(common_samples) > 0) {
        ancestry <- ancestry[common_samples, , drop = FALSE]
        majority_ancestry <- colnames(ancestry)[apply(ancestry, 1, which.max)]
        names(majority_ancestry) <- common_samples
        burden_df$ancestry_group <- majority_ancestry[burden_df$sample]
        cat(sprintf("  Ancestry: %d samples with ancestry assignments\n",
                    sum(!is.na(burden_df$ancestry_group))))
    }
}

# ============================================================================
# Statistical comparisons
# ============================================================================
cat("\n--- Statistical Comparisons ---\n")

comparison_results <- list()

# Helper function for Wilcoxon comparison
run_wilcoxon <- function(df, group_col, group1, group2, value_col, comparison_name) {
    g1_vals <- df[[value_col]][df[[group_col]] == group1]
    g2_vals <- df[[value_col]][df[[group_col]] == group2]

    g1_vals <- g1_vals[!is.na(g1_vals)]
    g2_vals <- g2_vals[!is.na(g2_vals)]

    if (length(g1_vals) < 3 || length(g2_vals) < 3) {
        cat(sprintf("  %s: insufficient samples (%s=%d, %s=%d), skipping\n",
                    comparison_name, group1, length(g1_vals), group2, length(g2_vals)))
        return(NULL)
    }

    test_res <- wilcox.test(g1_vals, g2_vals, exact = FALSE)

    result <- data.frame(
        comparison = comparison_name,
        value = value_col,
        group1 = group1,
        group2 = group2,
        n_group1 = length(g1_vals),
        n_group2 = length(g2_vals),
        median_group1 = median(g1_vals),
        median_group2 = median(g2_vals),
        wilcox_pvalue = test_res$p.value,
        stringsAsFactors = FALSE
    )
    cat(sprintf("  %s [%s]: %s median=%.1f (n=%d) vs %s median=%.1f (n=%d), p=%.4g\n",
                comparison_name, value_col, group1, result$median_group1,
                result$n_group1, group2, result$median_group2,
                result$n_group2, result$wilcox_pvalue))
    return(result)
}

# -----------------------------------------------------------------
# 1. Relapse vs non-relapse
# -----------------------------------------------------------------
relapse_cols <- c("relapse", "relapsed", "relapse_status", "event")
for (rcol in relapse_cols) {
    if (rcol %in% colnames(burden_df)) {
        cat(sprintf("\n  Comparing by %s...\n", rcol))
        groups <- unique(burden_df[[rcol]][!is.na(burden_df[[rcol]])])

        # Try to identify relapse/non-relapse groups
        relapse_group <- groups[grepl("yes|relapse|1|TRUE|recurrence", groups,
                                      ignore.case = TRUE)]
        no_relapse_group <- groups[grepl("no|non|0|FALSE|event.free", groups,
                                         ignore.case = TRUE)]

        if (length(relapse_group) >= 1 && length(no_relapse_group) >= 1) {
            for (val_col in c("total_burden", source_cols)) {
                result <- run_wilcoxon(burden_df, rcol,
                                       relapse_group[1], no_relapse_group[1],
                                       val_col, paste0("relapse_", val_col))
                if (!is.null(result)) {
                    comparison_results[[paste0("relapse_", val_col)]] <- result
                }
            }
        }
        break
    }
}

# -----------------------------------------------------------------
# 2. MRD-positive vs MRD-negative
# -----------------------------------------------------------------
mrd_cols <- c("MRD", "mrd", "mrd_status", "MRD_status")
for (mcol in mrd_cols) {
    if (mcol %in% colnames(burden_df)) {
        cat(sprintf("\n  Comparing by %s...\n", mcol))
        groups <- unique(burden_df[[mcol]][!is.na(burden_df[[mcol]])])

        pos_group <- groups[grepl("pos|\\+|1|TRUE|detected", groups,
                                   ignore.case = TRUE)]
        neg_group <- groups[grepl("neg|\\-|0|FALSE|undetected|not.detected", groups,
                                   ignore.case = TRUE)]

        if (length(pos_group) >= 1 && length(neg_group) >= 1) {
            for (val_col in c("total_burden", source_cols)) {
                result <- run_wilcoxon(burden_df, mcol,
                                       pos_group[1], neg_group[1],
                                       val_col, paste0("mrd_", val_col))
                if (!is.null(result)) {
                    comparison_results[[paste0("mrd_", val_col)]] <- result
                }
            }
        }
        break
    }
}

# -----------------------------------------------------------------
# 3. Ancestry group comparisons (Kruskal-Wallis)
# -----------------------------------------------------------------
if (!is.null(majority_ancestry) && "ancestry_group" %in% colnames(burden_df)) {
    cat("\n  Comparing by ancestry groups (Kruskal-Wallis)...\n")

    anc_valid <- burden_df %>%
        filter(!is.na(ancestry_group))

    anc_table <- table(anc_valid$ancestry_group)
    valid_groups <- names(anc_table[anc_table >= 3])

    if (length(valid_groups) >= 2) {
        anc_subset <- anc_valid %>% filter(ancestry_group %in% valid_groups)

        for (val_col in c("total_burden", source_cols)) {
            vals <- anc_subset[[val_col]]
            grps <- factor(anc_subset$ancestry_group)

            if (sum(!is.na(vals)) >= length(valid_groups) * 3) {
                kw_res <- kruskal.test(vals ~ grps)

                # Per-group medians
                group_medians <- anc_subset %>%
                    group_by(ancestry_group) %>%
                    summarise(
                        n = n(),
                        median_burden = median(!!sym(val_col), na.rm = TRUE),
                        .groups = "drop"
                    )

                kw_result <- data.frame(
                    comparison = paste0("ancestry_", val_col),
                    value = val_col,
                    group1 = paste(valid_groups, collapse = "/"),
                    group2 = "multi-group",
                    n_group1 = nrow(anc_subset),
                    n_group2 = length(valid_groups),
                    median_group1 = median(vals, na.rm = TRUE),
                    median_group2 = NA_real_,
                    wilcox_pvalue = kw_res$p.value,
                    stringsAsFactors = FALSE
                )

                comparison_results[[paste0("ancestry_", val_col)]] <- kw_result

                cat(sprintf("  Ancestry [%s]: KW p=%.4g across %d groups (%s)\n",
                            val_col, kw_res$p.value, length(valid_groups),
                            paste(valid_groups, collapse = ", ")))

                # Save per-group medians
                write.table(group_medians,
                            file.path(comp_dir,
                                      paste0("ancestry_", val_col, "_group_medians.tsv")),
                            sep = "\t", quote = FALSE, row.names = FALSE)
            }
        }
    } else {
        cat("  Fewer than 2 ancestry groups with sufficient samples, skipping\n")
    }
}

# Save all comparison results
if (length(comparison_results) > 0) {
    all_comparisons <- bind_rows(comparison_results)
    all_comparisons$padj <- p.adjust(all_comparisons$wilcox_pvalue, method = "BH")
    write.table(all_comparisons,
                file.path(comp_dir, "neoantigen_burden_comparisons.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("\n  Saved %d comparisons\n", nrow(all_comparisons)))
}

# ============================================================================
# Generate Plots
# ============================================================================
cat("\n--- Generating Plots ---\n")

# -----------------------------------------------------------------
# Boxplots for each comparison
# -----------------------------------------------------------------

# Relapse boxplot
for (rcol in relapse_cols) {
    if (rcol %in% colnames(burden_df)) {
        tryCatch({
            plot_df <- burden_df %>% filter(!is.na(!!sym(rcol)))
            if (nrow(plot_df) > 0) {
                p <- ggplot(plot_df,
                            aes(x = !!sym(rcol), y = total_burden, fill = !!sym(rcol))) +
                    geom_boxplot(outlier.size = 1) +
                    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
                    theme_minimal() +
                    labs(
                        title = "Neoantigen Burden: Relapse vs Non-Relapse",
                        x = "Relapse Status",
                        y = "Total Neoantigen Burden",
                        fill = "Group"
                    ) +
                    theme(legend.position = "none")

                ggsave(file.path(plot_dir, "burden_relapse_boxplot.pdf"),
                       p, width = 6, height = 5)
                cat("  Saved relapse boxplot\n")
            }
        }, error = function(e) {
            cat(sprintf("  WARNING: Relapse boxplot failed: %s\n", conditionMessage(e)))
        })
        break
    }
}

# MRD boxplot
for (mcol in mrd_cols) {
    if (mcol %in% colnames(burden_df)) {
        tryCatch({
            plot_df <- burden_df %>% filter(!is.na(!!sym(mcol)))
            if (nrow(plot_df) > 0) {
                p <- ggplot(plot_df,
                            aes(x = !!sym(mcol), y = total_burden, fill = !!sym(mcol))) +
                    geom_boxplot(outlier.size = 1) +
                    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
                    theme_minimal() +
                    labs(
                        title = "Neoantigen Burden: MRD-Positive vs MRD-Negative",
                        x = "MRD Status",
                        y = "Total Neoantigen Burden",
                        fill = "Group"
                    ) +
                    theme(legend.position = "none")

                ggsave(file.path(plot_dir, "burden_mrd_boxplot.pdf"),
                       p, width = 6, height = 5)
                cat("  Saved MRD boxplot\n")
            }
        }, error = function(e) {
            cat(sprintf("  WARNING: MRD boxplot failed: %s\n", conditionMessage(e)))
        })
        break
    }
}

# Ancestry boxplot
if ("ancestry_group" %in% colnames(burden_df)) {
    tryCatch({
        plot_df <- burden_df %>% filter(!is.na(ancestry_group))
        if (nrow(plot_df) > 0) {
            p <- ggplot(plot_df,
                        aes(x = ancestry_group, y = total_burden,
                            fill = ancestry_group)) +
                geom_boxplot(outlier.size = 1) +
                geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
                theme_minimal() +
                labs(
                    title = "Neoantigen Burden by Ancestry Group",
                    x = "Ancestry Group",
                    y = "Total Neoantigen Burden",
                    fill = "Ancestry"
                ) +
                theme(legend.position = "none",
                      axis.text.x = element_text(angle = 45, hjust = 1))

            ggsave(file.path(plot_dir, "burden_ancestry_boxplot.pdf"),
                   p, width = 8, height = 5)
            cat("  Saved ancestry boxplot\n")
        }
    }, error = function(e) {
        cat(sprintf("  WARNING: Ancestry boxplot failed: %s\n", conditionMessage(e)))
    })
}

# -----------------------------------------------------------------
# Stacked bar chart by neoantigen source
# -----------------------------------------------------------------
if (length(source_cols) > 1) {
    tryCatch({
        bar_data <- burden_df %>%
            select(sample, all_of(source_cols)) %>%
            pivot_longer(cols = all_of(source_cols),
                         names_to = "source", values_to = "count") %>%
            mutate(source = factor(source, levels = source_cols))

        # Order samples by total burden
        sample_order <- burden_df %>%
            arrange(desc(total_burden)) %>%
            pull(sample)
        bar_data$sample <- factor(bar_data$sample, levels = sample_order)

        p_stack <- ggplot(bar_data, aes(x = sample, y = count, fill = source)) +
            geom_bar(stat = "identity") +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                legend.position = "bottom"
            ) +
            scale_fill_brewer(palette = "Set2") +
            labs(
                title = "Neoantigen Burden by Source",
                x = "Sample",
                y = "Neoantigen Count",
                fill = "Source"
            )

        # Adjust width based on number of samples
        plot_width <- max(8, nrow(burden_df) * 0.3)
        ggsave(file.path(plot_dir, "burden_stacked_barchart.pdf"),
               p_stack, width = min(plot_width, 24), height = 6)
        cat("  Saved stacked bar chart\n")
    }, error = function(e) {
        cat(sprintf("  WARNING: Stacked bar chart failed: %s\n", conditionMessage(e)))
    })
}

# ============================================================================
# NESTED STRATIFICATION: MRD-negative relapse vs non-relapse + ancestry
# ============================================================================
cat("\n--- Nested Stratification: MRD-Negative Strata ---\n")

# Identify MRD and relapse columns in the data
mrd_col_found <- NULL
for (mcol in mrd_cols) {
    if (mcol %in% colnames(burden_df)) { mrd_col_found <- mcol; break }
}
relapse_col_found <- NULL
for (rcol in relapse_cols) {
    if (rcol %in% colnames(burden_df)) { relapse_col_found <- rcol; break }
}

if (!is.null(mrd_col_found) && !is.null(relapse_col_found)) {
    mrd_groups <- unique(burden_df[[mrd_col_found]][!is.na(burden_df[[mrd_col_found]])])
    neg_group_val <- mrd_groups[grepl("neg|\\-|0|FALSE|undetected|not.detected",
                                      mrd_groups, ignore.case = TRUE)]

    if (length(neg_group_val) >= 1) {
        mrd_neg <- burden_df %>% filter(!!sym(mrd_col_found) == neg_group_val[1],
                                         !is.na(!!sym(relapse_col_found)))
        cat(sprintf("  MRD-negative subset: %d samples\n", nrow(mrd_neg)))

        relapse_groups <- unique(mrd_neg[[relapse_col_found]])
        relapse_grp <- relapse_groups[grepl("yes|relapse|1|TRUE|recurrence",
                                             relapse_groups, ignore.case = TRUE)]
        no_relapse_grp <- relapse_groups[grepl("no|non|0|FALSE|event.free",
                                                relapse_groups, ignore.case = TRUE)]

        if (nrow(mrd_neg) >= 4 && length(relapse_grp) >= 1 && length(no_relapse_grp) >= 1) {
            nested_dir <- file.path(comp_dir, "mrd_neg_stratified")
            nested_plot_dir <- file.path(plot_dir, "mrd_neg_stratified")
            dir.create(nested_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(nested_plot_dir, recursive = TRUE, showWarnings = FALSE)

            # --- Relapse vs non-relapse within MRD-negative ---
            nested_results <- list()
            for (val_col in c("total_burden", source_cols)) {
                result <- run_wilcoxon(mrd_neg, relapse_col_found,
                                       relapse_grp[1], no_relapse_grp[1],
                                       val_col, paste0("mrd_neg_relapse_", val_col))
                if (!is.null(result)) {
                    nested_results[[paste0("mrd_neg_relapse_", val_col)]] <- result
                }
            }
            if (length(nested_results) > 0) {
                nested_comp <- bind_rows(nested_results)
                nested_comp$padj <- p.adjust(nested_comp$wilcox_pvalue, method = "BH")
                write.table(nested_comp,
                            file.path(nested_dir, "relapse_within_mrd_neg.tsv"),
                            sep = "\t", quote = FALSE, row.names = FALSE)
            }

            # Boxplot: total burden by relapse within MRD-neg
            tryCatch({
                p_nested <- ggplot(mrd_neg,
                                   aes(x = !!sym(relapse_col_found), y = total_burden,
                                       fill = !!sym(relapse_col_found))) +
                    geom_boxplot(outlier.shape = NA) +
                    geom_jitter(width = 0.2, alpha = 0.5) +
                    labs(title = "MRD-Negative: Neoantigen Burden by Relapse",
                         x = "", y = "Total Neoantigens") +
                    theme_bw() + theme(legend.position = "none")
                ggsave(file.path(nested_plot_dir, "mrd_neg_burden_by_relapse.pdf"),
                       p_nested, width = 6, height = 5)
            }, error = function(e) cat(sprintf("  WARNING: nested boxplot failed: %s\n",
                                                conditionMessage(e))))

            # Faceted by neoantigen type
            if (length(source_cols) > 1) {
                tryCatch({
                    mrd_neg_long <- mrd_neg %>%
                        select(sample, !!sym(relapse_col_found), all_of(source_cols)) %>%
                        pivot_longer(cols = all_of(source_cols),
                                     names_to = "neoantigen_type", values_to = "count")

                    p_types <- ggplot(mrd_neg_long,
                                      aes(x = !!sym(relapse_col_found), y = count,
                                          fill = !!sym(relapse_col_found))) +
                        geom_boxplot(outlier.shape = NA) +
                        geom_jitter(width = 0.2, alpha = 0.5) +
                        facet_wrap(~neoantigen_type, scales = "free_y") +
                        labs(title = "MRD-Negative: Neoantigen Types by Relapse",
                             x = "", y = "Count") +
                        theme_bw() + theme(legend.position = "none")
                    ggsave(file.path(nested_plot_dir, "mrd_neg_types_by_relapse.pdf"),
                           p_types, width = 10, height = 8)
                }, error = function(e) cat(sprintf("  WARNING: type facet plot failed: %s\n",
                                                    conditionMessage(e))))
            }

            # --- Ancestry category within MRD-negative relapse strata ---
            if ("ancestry_group" %in% colnames(mrd_neg)) {
                for (rel_val in c(relapse_grp[1], no_relapse_grp[1])) {
                    sub <- mrd_neg %>%
                        filter(!!sym(relapse_col_found) == rel_val, !is.na(ancestry_group))
                    valid_anc <- names(which(table(sub$ancestry_group) >= 3))

                    if (nrow(sub) >= 4 && length(valid_anc) >= 2) {
                        sub_filt <- sub %>% filter(ancestry_group %in% valid_anc)
                        anc_nested <- list()
                        for (val_col in c("total_burden", source_cols)) {
                            kt <- tryCatch(
                                kruskal.test(as.formula(paste(val_col, "~ ancestry_group")),
                                             data = sub_filt),
                                error = function(e) NULL)
                            if (!is.null(kt)) {
                                grp_medians <- sub_filt %>%
                                    group_by(ancestry_group) %>%
                                    summarise(n = n(), median = median(!!sym(val_col), na.rm = TRUE),
                                              .groups = "drop")
                                anc_nested[[val_col]] <- data.frame(
                                    stratum = paste0("MRD_neg_", rel_val),
                                    value = val_col,
                                    kw_pvalue = kt$p.value,
                                    n_groups = length(valid_anc),
                                    stringsAsFactors = FALSE
                                )
                                cat(sprintf("  MRD-neg/%s ancestry [%s]: KW p=%.4g\n",
                                            rel_val, val_col, kt$p.value))
                            }
                        }
                        if (length(anc_nested) > 0) {
                            anc_df <- bind_rows(anc_nested)
                            anc_df$padj <- p.adjust(anc_df$kw_pvalue, method = "BH")
                            write.table(anc_df,
                                        file.path(nested_dir,
                                                  paste0("ancestry_in_mrd_neg_", rel_val, ".tsv")),
                                        sep = "\t", quote = FALSE, row.names = FALSE)
                        }
                    }
                }

                # Ancestry x burden within MRD-neg, faceted by relapse
                tryCatch({
                    p_anc_rel <- ggplot(mrd_neg %>% filter(!is.na(ancestry_group)),
                                        aes(x = ancestry_group, y = total_burden,
                                            fill = ancestry_group)) +
                        geom_boxplot(outlier.shape = NA) +
                        geom_jitter(width = 0.2, alpha = 0.5) +
                        facet_wrap(as.formula(paste("~", relapse_col_found))) +
                        labs(title = "MRD-Neg: Burden by Ancestry & Relapse",
                             x = "", y = "Total Neoantigens") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                              legend.position = "none")
                    ggsave(file.path(nested_plot_dir, "mrd_neg_ancestry_by_relapse.pdf"),
                           p_anc_rel, width = 10, height = 5)
                }, error = function(e) cat(sprintf("  WARNING: ancestry x relapse plot failed: %s\n",
                                                    conditionMessage(e))))
            }

            # --- Ancestry proportion correlations within MRD-negative strata ---
            if (!is.null(ancestry)) {
                anc_prop_cols <- colnames(ancestry)
                cor_results <- list()
                for (rel_val in c(relapse_grp[1], no_relapse_grp[1])) {
                    sub <- mrd_neg %>% filter(!!sym(relapse_col_found) == rel_val)
                    sub_samples <- intersect(sub$sample, rownames(ancestry))
                    if (length(sub_samples) >= 5) {
                        for (acol in anc_prop_cols) {
                            anc_vals <- ancestry[sub_samples, acol]
                            for (bcol in c("total_burden", source_cols)) {
                                burd_vals <- sub[[bcol]][match(sub_samples, sub$sample)]
                                ct <- tryCatch(cor.test(anc_vals, burd_vals,
                                                        method = "spearman", exact = FALSE),
                                               error = function(e) NULL)
                                if (!is.null(ct)) {
                                    cor_results[[length(cor_results) + 1]] <- data.frame(
                                        stratum = paste0("MRD_neg_", rel_val),
                                        ancestry_variable = acol,
                                        neoantigen_type = bcol,
                                        spearman_rho = ct$estimate,
                                        p_value = ct$p.value,
                                        n = length(sub_samples),
                                        stringsAsFactors = FALSE
                                    )
                                }
                            }
                        }
                    }
                }
                if (length(cor_results) > 0) {
                    cor_df <- bind_rows(cor_results)
                    cor_df$padj <- p.adjust(cor_df$p_value, method = "BH")
                    cor_df <- cor_df %>% arrange(padj)
                    write.table(cor_df,
                                file.path(nested_dir, "ancestry_proportion_correlations.tsv"),
                                sep = "\t", quote = FALSE, row.names = FALSE)
                    cat(sprintf("  Saved %d ancestry proportion correlations\n", nrow(cor_df)))

                    # Scatter for top correlation
                    best <- cor_df[1, ]
                    if (!is.na(best$padj)) {
                        rel_grp_best <- gsub("MRD_neg_", "", best$stratum)
                        sub_best <- mrd_neg %>% filter(!!sym(relapse_col_found) == rel_grp_best)
                        sub_samples <- intersect(sub_best$sample, rownames(ancestry))
                        plot_df <- data.frame(
                            anc = ancestry[sub_samples, best$ancestry_variable],
                            burden = sub_best[[best$neoantigen_type]][match(sub_samples, sub_best$sample)]
                        )
                        tryCatch({
                            p_cor <- ggplot(plot_df, aes(x = anc, y = burden)) +
                                geom_point(alpha = 0.7) +
                                geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
                                labs(title = sprintf("MRD-Neg/%s: %s vs %s (rho=%.2f, p=%.2e)",
                                                     rel_grp_best, best$ancestry_variable,
                                                     best$neoantigen_type,
                                                     best$spearman_rho, best$p_value),
                                     x = best$ancestry_variable, y = best$neoantigen_type) +
                                theme_bw()
                            ggsave(file.path(nested_plot_dir, "top_ancestry_correlation.pdf"),
                                   p_cor, width = 7, height = 5)
                        }, error = function(e) cat(sprintf("  WARNING: correlation plot failed: %s\n",
                                                            conditionMessage(e))))
                    }
                }
            }
        }
    }
} else {
    cat("  MRD and/or relapse columns not found, skipping nested stratification\n")
}

# ============================================================================
# Save annotated burden table
# ============================================================================
cat("\n--- Saving Annotated Burden Table ---\n")

# Select relevant columns for output
output_cols <- c("sample", source_cols, "total_burden")
if ("ancestry_group" %in% colnames(burden_df)) {
    output_cols <- c(output_cols, "ancestry_group")
}
for (rcol in relapse_cols) {
    if (rcol %in% colnames(burden_df)) {
        output_cols <- c(output_cols, rcol)
        break
    }
}
for (mcol in mrd_cols) {
    if (mcol %in% colnames(burden_df)) {
        output_cols <- c(output_cols, mcol)
        break
    }
}

output_cols <- intersect(output_cols, colnames(burden_df))
burden_output <- burden_df[, output_cols, drop = FALSE] %>%
    arrange(desc(total_burden))

write.table(burden_output,
            file.path(opt$output_dir, "annotated_neoantigen_burden.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("  Saved annotated burden table: %d samples x %d columns\n",
            nrow(burden_output), ncol(burden_output)))

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Neoantigen Burden Analysis Complete ===\n")
cat(sprintf("Output directory: %s\n", opt$output_dir))
cat("Key outputs:\n")
cat(sprintf("  - Annotated burden table: %s\n",
            file.path(opt$output_dir, "annotated_neoantigen_burden.tsv")))
cat(sprintf("  - Statistical comparisons: %s\n", comp_dir))
cat(sprintf("  - Plots: %s\n", plot_dir))
cat(sprintf("  Median total burden: %.1f\n", median(burden_df$total_burden, na.rm = TRUE)))
cat(sprintf("  Range: %d - %d neoantigens\n",
            min(burden_df$total_burden, na.rm = TRUE),
            max(burden_df$total_burden, na.rm = TRUE)))
