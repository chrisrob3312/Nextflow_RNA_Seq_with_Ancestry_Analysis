#!/usr/bin/env Rscript

# ============================================================================
# Bisbee Diff: Beta-binomial differential splicing analysis
#
# Reads Bisbee prep output TSVs (one per event type) and fits a beta-binomial
# GLM (VGAM package) to test for differential splicing between two conditions.
# Produces a combined results table and a volcano plot.
#
# Input:  Per-event-type TSVs from bisbee_prep.py with columns:
#         event_id, gene, sample_id, inclusion_count, exclusion_count, psi
# Output: differential_splicing.tsv, volcano_plot.pdf
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(VGAM)
    library(ggplot2)
    library(parallel)
})

# ---- Parse arguments --------------------------------------------------------

option_list <- list(
    make_option("--input-dir", type = "character",
                help = "Directory containing Bisbee prep output TSVs"),
    make_option("--conditions", type = "character",
                help = "TSV file mapping sample_id to condition (columns: sample_id, condition)"),
    make_option("--contrast", type = "character",
                help = "Contrast specification as 'group1_vs_group2' (e.g., 'tumor_vs_normal')"),
    make_option("--output-dir", type = "character", default = "bisbee_diff_results",
                help = "Output directory for results [default: %default]"),
    make_option("--threads", type = "integer", default = 1,
                help = "Number of parallel threads [default: %default]"),
    make_option("--min-samples", type = "integer", default = 3,
                help = "Minimum samples per group with non-zero counts [default: %default]"),
    make_option("--min-total-count", type = "integer", default = 10,
                help = "Minimum total count (inc + exc) per event-sample [default: %default]"),
    make_option("--padj-threshold", type = "double", default = 0.05,
                help = "Adjusted p-value threshold for significance [default: %default]"),
    make_option("--dpsi-threshold", type = "double", default = 0.1,
                help = "Minimum absolute delta-PSI for significance [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate inputs --------------------------------------------------------

if (is.null(opt$`input-dir`) || !dir.exists(opt$`input-dir`)) {
    stop("ERROR: --input-dir must be an existing directory containing Bisbee prep TSVs")
}

if (is.null(opt$conditions) || !file.exists(opt$conditions)) {
    stop("ERROR: --conditions must point to an existing TSV file")
}

if (is.null(opt$contrast)) {
    stop("ERROR: --contrast must be specified as 'group1_vs_group2'")
}

# Parse contrast
contrast_parts <- strsplit(opt$contrast, "_vs_")[[1]]
if (length(contrast_parts) != 2) {
    stop("ERROR: --contrast must be formatted as 'group1_vs_group2'")
}
group1_name <- trimws(contrast_parts[1])
group2_name <- trimws(contrast_parts[2])

cat("Contrast:", group1_name, "vs", group2_name, "\n")

# Create output directory
dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- Load condition mapping -------------------------------------------------

conditions <- read.delim(opt$conditions, stringsAsFactors = FALSE)
if (!all(c("sample_id", "condition") %in% colnames(conditions))) {
    stop("ERROR: Conditions file must have columns 'sample_id' and 'condition'")
}

group1_samples <- conditions$sample_id[conditions$condition == group1_name]
group2_samples <- conditions$sample_id[conditions$condition == group2_name]

if (length(group1_samples) == 0) {
    stop(paste("ERROR: No samples found for group:", group1_name))
}
if (length(group2_samples) == 0) {
    stop(paste("ERROR: No samples found for group:", group2_name))
}

cat("Group 1 (", group1_name, "):", length(group1_samples), "samples\n")
cat("Group 2 (", group2_name, "):", length(group2_samples), "samples\n")

# ---- Load all event-type TSVs -----------------------------------------------

tsv_files <- list.files(opt$`input-dir`, pattern = "_counts\\.tsv$", full.names = TRUE)

if (length(tsv_files) == 0) {
    cat("WARNING: No *_counts.tsv files found in", opt$`input-dir`, "\n")
    # Write empty output
    empty_df <- data.frame(
        event_id = character(), gene = character(), event_type = character(),
        deltaPSI = numeric(), pvalue = numeric(), padj = numeric(),
        mean_psi_group1 = numeric(), mean_psi_group2 = numeric(),
        stringsAsFactors = FALSE
    )
    write.table(empty_df, file.path(opt$`output-dir`, "differential_splicing.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("No input data. Empty output written.\n")
    quit(save = "no", status = 0)
}

cat("Loading", length(tsv_files), "event-type count file(s)...\n")

all_data <- do.call(rbind, lapply(tsv_files, function(f) {
    # Derive event type from filename (e.g., exon_skip_counts.tsv -> exon_skip)
    event_type <- sub("_counts\\.tsv$", "", basename(f))
    df <- tryCatch(
        read.delim(f, stringsAsFactors = FALSE),
        error = function(e) {
            cat("WARNING: Could not read", f, ":", e$message, "\n")
            return(NULL)
        }
    )
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$event_type <- event_type
    return(df)
}))

if (is.null(all_data) || nrow(all_data) == 0) {
    cat("WARNING: All input files are empty.\n")
    empty_df <- data.frame(
        event_id = character(), gene = character(), event_type = character(),
        deltaPSI = numeric(), pvalue = numeric(), padj = numeric(),
        mean_psi_group1 = numeric(), mean_psi_group2 = numeric(),
        stringsAsFactors = FALSE
    )
    write.table(empty_df, file.path(opt$`output-dir`, "differential_splicing.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    quit(save = "no", status = 0)
}

# Filter to samples present in conditions
all_samples <- union(group1_samples, group2_samples)
all_data <- all_data[all_data$sample_id %in% all_samples, ]

# Convert PSI to numeric (handle "NA" strings)
all_data$psi <- suppressWarnings(as.numeric(all_data$psi))

cat("Total rows after sample filtering:", nrow(all_data), "\n")
cat("Unique events:", length(unique(all_data$event_id)), "\n")

# ---- Beta-binomial differential splicing test per event ----------------------

fit_betabinom_event <- function(event_df, group1_samples, group2_samples,
                                 min_samples, min_total_count) {
    #' Fit a beta-binomial GLM for a single splicing event.
    #'
    #' @param event_df  Data frame with rows for one event across samples.
    #' @param group1_samples  Character vector of sample IDs in group 1.
    #' @param group2_samples  Character vector of sample IDs in group 2.
    #' @return Named list with deltaPSI, pvalue, mean_psi_group1, mean_psi_group2,
    #'         or NULL if the event cannot be tested.

    event_id <- event_df$event_id[1]
    gene <- event_df$gene[1]
    event_type <- event_df$event_type[1]

    # Assign group labels
    event_df$group <- NA_character_
    event_df$group[event_df$sample_id %in% group1_samples] <- "group1"
    event_df$group[event_df$sample_id %in% group2_samples] <- "group2"
    event_df <- event_df[!is.na(event_df$group), ]

    # Require minimum total count
    event_df$total <- event_df$inclusion_count + event_df$exclusion_count
    event_df <- event_df[event_df$total >= min_total_count, ]

    # Check minimum samples per group
    n_g1 <- sum(event_df$group == "group1")
    n_g2 <- sum(event_df$group == "group2")

    if (n_g1 < min_samples || n_g2 < min_samples) {
        return(NULL)
    }

    # Compute mean PSI per group
    g1_psi <- event_df$psi[event_df$group == "group1"]
    g2_psi <- event_df$psi[event_df$group == "group2"]
    # Remove NA PSI values for mean computation
    mean_psi_g1 <- mean(g1_psi[!is.na(g1_psi)])
    mean_psi_g2 <- mean(g2_psi[!is.na(g2_psi)])
    delta_psi <- mean_psi_g1 - mean_psi_g2

    # Fit beta-binomial GLM: inclusion_count modeled as beta-binomial with
    # total = inclusion + exclusion, testing group effect
    event_df$group <- factor(event_df$group, levels = c("group2", "group1"))

    pvalue <- NA_real_

    tryCatch({
        # Full model with group effect
        fit_full <- vglm(
            cbind(inclusion_count, exclusion_count) ~ group,
            family = betabinomial,
            data = event_df
        )

        # Null model (intercept only)
        fit_null <- vglm(
            cbind(inclusion_count, exclusion_count) ~ 1,
            family = betabinomial,
            data = event_df
        )

        # Likelihood ratio test
        lr_stat <- 2 * (logLik(fit_full) - logLik(fit_null))
        df_diff <- length(coef(fit_full)) - length(coef(fit_null))

        if (df_diff > 0 && is.finite(lr_stat) && lr_stat >= 0) {
            pvalue <- pchisq(as.numeric(lr_stat), df = df_diff, lower.tail = FALSE)
        }
    }, error = function(e) {
        # Model fitting failed; leave pvalue as NA
    }, warning = function(w) {
        # Suppress convergence warnings but continue
        invokeRestart("muffleWarning")
    })

    return(list(
        event_id = event_id,
        gene = gene,
        event_type = event_type,
        deltaPSI = delta_psi,
        pvalue = pvalue,
        mean_psi_group1 = mean_psi_g1,
        mean_psi_group2 = mean_psi_g2
    ))
}

# Split data by event
event_list <- split(all_data, all_data$event_id)
n_events <- length(event_list)
cat("Testing", n_events, "events for differential splicing...\n")

# Run beta-binomial tests (optionally in parallel)
if (opt$threads > 1 && .Platform$OS.type == "unix") {
    cat("Using", opt$threads, "threads\n")
    results_list <- mclapply(
        event_list,
        fit_betabinom_event,
        group1_samples = group1_samples,
        group2_samples = group2_samples,
        min_samples = opt$`min-samples`,
        min_total_count = opt$`min-total-count`,
        mc.cores = opt$threads
    )
} else {
    results_list <- lapply(
        event_list,
        fit_betabinom_event,
        group1_samples = group1_samples,
        group2_samples = group2_samples,
        min_samples = opt$`min-samples`,
        min_total_count = opt$`min-total-count`
    )
}

# Remove NULL entries (events that could not be tested)
results_list <- results_list[!sapply(results_list, is.null)]

if (length(results_list) == 0) {
    cat("WARNING: No events passed filtering criteria.\n")
    empty_df <- data.frame(
        event_id = character(), gene = character(), event_type = character(),
        deltaPSI = numeric(), pvalue = numeric(), padj = numeric(),
        mean_psi_group1 = numeric(), mean_psi_group2 = numeric(),
        stringsAsFactors = FALSE
    )
    write.table(empty_df, file.path(opt$`output-dir`, "differential_splicing.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    quit(save = "no", status = 0)
}

# Combine into data frame
results_df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(x, stringsAsFactors = FALSE)
}))

# Multiple testing correction (BH)
results_df$padj <- p.adjust(results_df$pvalue, method = "BH")

# Sort by adjusted p-value
results_df <- results_df[order(results_df$padj, -abs(results_df$deltaPSI)), ]

# Round numeric columns for readability
results_df$deltaPSI <- round(results_df$deltaPSI, 6)
results_df$mean_psi_group1 <- round(results_df$mean_psi_group1, 6)
results_df$mean_psi_group2 <- round(results_df$mean_psi_group2, 6)

# Reorder columns
results_df <- results_df[, c("event_id", "gene", "event_type", "deltaPSI",
                              "pvalue", "padj", "mean_psi_group1", "mean_psi_group2")]

# ---- Write results ----------------------------------------------------------

output_tsv <- file.path(opt$`output-dir`, "differential_splicing.tsv")
write.table(results_df, output_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

n_sig <- sum(results_df$padj < opt$`padj-threshold` &
             abs(results_df$deltaPSI) >= opt$`dpsi-threshold`, na.rm = TRUE)
n_tested <- sum(!is.na(results_df$pvalue))

cat("\nResults summary:\n")
cat("  Events tested:", n_tested, "\n")
cat("  Significant (padj <", opt$`padj-threshold`,
    ", |deltaPSI| >=", opt$`dpsi-threshold`, "):", n_sig, "\n")
cat("  Output:", output_tsv, "\n")

# ---- Volcano plot ------------------------------------------------------------

tryCatch({
    plot_df <- results_df[!is.na(results_df$pvalue) & !is.na(results_df$padj), ]

    if (nrow(plot_df) > 0) {
        plot_df$neg_log10_padj <- -log10(plot_df$padj)
        # Cap extreme values for plotting
        plot_df$neg_log10_padj <- pmin(plot_df$neg_log10_padj, 50)

        # Classify significance
        plot_df$significance <- "NS"
        plot_df$significance[abs(plot_df$deltaPSI) >= opt$`dpsi-threshold` &
                             plot_df$padj < opt$`padj-threshold`] <- "Significant"
        plot_df$significance[abs(plot_df$deltaPSI) >= opt$`dpsi-threshold` &
                             plot_df$padj >= opt$`padj-threshold`] <- "deltaPSI only"
        plot_df$significance[abs(plot_df$deltaPSI) < opt$`dpsi-threshold` &
                             plot_df$padj < opt$`padj-threshold`] <- "p-value only"

        # Label top significant events
        top_events <- head(
            plot_df[plot_df$significance == "Significant", ],
            min(20, sum(plot_df$significance == "Significant"))
        )

        volcano_colors <- c(
            "NS" = "grey70",
            "deltaPSI only" = "steelblue",
            "p-value only" = "darkorange",
            "Significant" = "firebrick"
        )

        p <- ggplot(plot_df, aes(x = deltaPSI, y = neg_log10_padj, color = significance)) +
            geom_point(alpha = 0.6, size = 1.5) +
            scale_color_manual(values = volcano_colors) +
            geom_hline(yintercept = -log10(opt$`padj-threshold`),
                       linetype = "dashed", color = "grey40") +
            geom_vline(xintercept = c(-opt$`dpsi-threshold`, opt$`dpsi-threshold`),
                       linetype = "dashed", color = "grey40") +
            labs(
                title = paste("Differential Splicing:", opt$contrast),
                subtitle = paste(n_sig, "significant events out of", n_tested, "tested"),
                x = expression(Delta * PSI),
                y = expression(-log[10](p[adj])),
                color = "Significance"
            ) +
            theme_bw(base_size = 12) +
            theme(
                legend.position = "bottom",
                plot.title = element_text(face = "bold")
            )

        # Add gene labels for top events
        if (nrow(top_events) > 0) {
            p <- p + ggplot2::annotate(
                "text",
                x = top_events$deltaPSI,
                y = top_events$neg_log10_padj,
                label = top_events$gene,
                size = 2.5, fontface = "italic",
                hjust = -0.1, vjust = -0.3
            )
        }

        plot_path <- file.path(opt$`output-dir`, "volcano_plot.pdf")
        ggsave(plot_path, p, width = 10, height = 8)
        cat("Volcano plot saved:", plot_path, "\n")
    } else {
        cat("No non-NA results to plot.\n")
    }
}, error = function(e) {
    cat("WARNING: Volcano plot generation failed:", e$message, "\n")
})

cat("\nBisbee differential splicing analysis complete.\n")
