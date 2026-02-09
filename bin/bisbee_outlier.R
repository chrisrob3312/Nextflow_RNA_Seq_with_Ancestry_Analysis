#!/usr/bin/env Rscript

# ============================================================================
# Bisbee Outlier: Per-sample splicing outlier detection
#
# For each sample, identifies splicing events with abnormally high or low PSI
# using a beta-binomial model fitted to the population distribution. Events
# are scored with z-scores against the population and flagged as outliers.
#
# Input:  Per-event-type TSVs from bisbee_prep.py
# Output: outlier_events.tsv (per-event outlier calls)
#         outlier_summary.tsv (per-sample outlier counts)
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(VGAM)
})

# ---- Parse arguments --------------------------------------------------------

option_list <- list(
    make_option("--input-dir", type = "character",
                help = "Directory containing Bisbee prep output TSVs"),
    make_option("--output-dir", type = "character", default = "bisbee_outlier_results",
                help = "Output directory for outlier results [default: %default]"),
    make_option("--fdr-threshold", type = "double", default = 0.05,
                help = "FDR threshold for outlier significance [default: %default]"),
    make_option("--zscore-threshold", type = "double", default = 2.0,
                help = "Absolute z-score threshold for outlier detection [default: %default]"),
    make_option("--min-samples", type = "integer", default = 10,
                help = "Minimum number of samples with data to fit population model [default: %default]"),
    make_option("--min-total-count", type = "integer", default = 10,
                help = "Minimum total junction count per event-sample pair [default: %default]"),
    make_option("--threads", type = "integer", default = 1,
                help = "Number of parallel threads [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate inputs --------------------------------------------------------

if (is.null(opt$`input-dir`) || !dir.exists(opt$`input-dir`)) {
    stop("ERROR: --input-dir must be an existing directory containing Bisbee prep TSVs")
}

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- Load all event-type TSVs -----------------------------------------------

tsv_files <- list.files(opt$`input-dir`, pattern = "_counts\\.tsv$", full.names = TRUE)

if (length(tsv_files) == 0) {
    cat("WARNING: No *_counts.tsv files found in", opt$`input-dir`, "\n")
    # Write empty outputs
    write_empty_outputs(opt$`output-dir`)
    quit(save = "no", status = 0)
}

cat("Loading", length(tsv_files), "event-type count file(s)...\n")

all_data <- do.call(rbind, lapply(tsv_files, function(f) {
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
    write_empty_outputs(opt$`output-dir`)
    quit(save = "no", status = 0)
}

# Convert PSI to numeric
all_data$psi <- suppressWarnings(as.numeric(all_data$psi))

# Compute total counts
all_data$total_count <- all_data$inclusion_count + all_data$exclusion_count

# Filter by minimum total count
all_data <- all_data[all_data$total_count >= opt$`min-total-count`, ]

n_events <- length(unique(all_data$event_id))
n_samples <- length(unique(all_data$sample_id))
cat("After filtering: ", nrow(all_data), " observations (",
    n_events, " events x ", n_samples, " samples)\n", sep = "")

if (nrow(all_data) == 0) {
    cat("WARNING: No data remaining after count filtering.\n")
    write_empty_outputs(opt$`output-dir`)
    quit(save = "no", status = 0)
}

# ---- Fit population beta-binomial model per event ---------------------------

#' Fit a beta-binomial model to the population distribution of PSI for one event,
#' then score each sample as an outlier.
#'
#' @param event_df  Data frame with rows for one event across all samples.
#' @param min_samples Minimum number of samples required.
#' @param zscore_threshold Absolute z-score threshold.
#' @return Data frame with outlier scores per sample for this event, or NULL.
fit_event_outliers <- function(event_df, min_samples, zscore_threshold) {

    event_id <- event_df$event_id[1]
    gene <- event_df$gene[1]
    event_type <- event_df$event_type[1]

    # Need enough samples for reliable population estimate
    if (nrow(event_df) < min_samples) {
        return(NULL)
    }

    # Remove rows with NA PSI
    event_df <- event_df[!is.na(event_df$psi), ]
    if (nrow(event_df) < min_samples) {
        return(NULL)
    }

    # Fit beta-binomial to population using method of moments as fallback
    # and VGAM MLE as primary method
    mu_hat <- NA_real_
    rho_hat <- NA_real_
    fit_success <- FALSE

    tryCatch({
        # VGAM beta-binomial fit (intercept-only model)
        fit <- vglm(
            cbind(inclusion_count, exclusion_count) ~ 1,
            family = betabinomial,
            data = event_df
        )

        # Extract mu (mean) and rho (overdispersion) from VGAM parameterization
        # VGAM betabinomial uses logit link for mu and log link for rho
        coefs <- Coef(fit)
        mu_hat <- coefs[1]   # Mean PSI (on probability scale after inverse-logit)
        rho_hat <- coefs[2]  # Overdispersion parameter

        # Validate parameters
        if (is.finite(mu_hat) && is.finite(rho_hat) && mu_hat > 0 && mu_hat < 1 && rho_hat > 0) {
            fit_success <- TRUE
        }
    }, error = function(e) {
        # MLE failed; fall through to method of moments
    }, warning = function(w) {
        invokeRestart("muffleWarning")
    })

    if (!fit_success) {
        # Method of moments fallback
        psi_vals <- event_df$psi
        psi_mean <- mean(psi_vals, na.rm = TRUE)
        psi_var <- var(psi_vals, na.rm = TRUE)

        # Avoid edge cases
        if (is.na(psi_mean) || is.na(psi_var) || psi_var == 0 ||
            psi_mean <= 0 || psi_mean >= 1) {
            return(NULL)
        }

        mu_hat <- psi_mean
        # Estimate rho from variance: Var(PSI) ~ mu*(1-mu)*(1+rho)/(1+n*rho) for beta-binomial
        # Simplified: rho ~ (Var/mu/(1-mu) - 1/n) / (1 - 1/n) where n is mean total count
        n_mean <- mean(event_df$total_count)
        theoretical_var <- mu_hat * (1 - mu_hat) / n_mean
        if (psi_var > theoretical_var) {
            rho_hat <- min((psi_var / (mu_hat * (1 - mu_hat)) - 1 / n_mean) /
                          (1 - 1 / n_mean), 0.99)
            rho_hat <- max(rho_hat, 0.001)
        } else {
            rho_hat <- 0.001
        }
    }

    # Compute z-scores for each sample
    # Under beta-binomial: E[PSI] = mu, Var[PSI] approx mu*(1-mu)*(1+(n-1)*rho)/n
    # where n is the total count for that sample
    results <- data.frame(
        event_id = character(),
        gene = character(),
        event_type = character(),
        sample_id = character(),
        psi = numeric(),
        population_mean = numeric(),
        zscore = numeric(),
        pvalue = numeric(),
        direction = character(),
        stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(event_df))) {
        row <- event_df[i, ]
        n_i <- row$total_count
        psi_i <- row$psi

        if (is.na(psi_i) || n_i < 1) next

        # Variance of PSI under beta-binomial
        var_psi <- mu_hat * (1 - mu_hat) * (1 + (n_i - 1) * rho_hat) / n_i
        if (var_psi <= 0 || !is.finite(var_psi)) next

        sd_psi <- sqrt(var_psi)
        zscore <- (psi_i - mu_hat) / sd_psi

        if (!is.finite(zscore)) next

        # Two-sided p-value from standard normal
        pval <- 2 * pnorm(-abs(zscore))

        # Only report if above z-score threshold (pre-filter before FDR)
        if (abs(zscore) >= zscore_threshold) {
            direction <- ifelse(zscore > 0, "high_PSI", "low_PSI")

            results <- rbind(results, data.frame(
                event_id = event_id,
                gene = gene,
                event_type = event_type,
                sample_id = row$sample_id,
                psi = round(psi_i, 6),
                population_mean = round(mu_hat, 6),
                zscore = round(zscore, 4),
                pvalue = pval,
                direction = direction,
                stringsAsFactors = FALSE
            ))
        }
    }

    if (nrow(results) == 0) return(NULL)

    return(results)
}


#' Write empty output files when there is no data to process.
write_empty_outputs <- function(output_dir) {
    # Empty outlier events
    events_header <- c("event_id", "gene", "event_type", "sample_id", "psi",
                       "population_mean", "zscore", "pvalue", "padj", "direction")
    write.table(
        data.frame(matrix(ncol = length(events_header), nrow = 0,
                          dimnames = list(NULL, events_header))),
        file.path(output_dir, "outlier_events.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    # Empty outlier summary
    summary_header <- c("sample_id", "n_outlier_events", "n_high_psi", "n_low_psi",
                        "n_event_types_affected", "top_outlier_gene")
    write.table(
        data.frame(matrix(ncol = length(summary_header), nrow = 0,
                          dimnames = list(NULL, summary_header))),
        file.path(output_dir, "outlier_summary.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    cat("Empty output files written to", output_dir, "\n")
}


# ---- Run outlier detection per event ----------------------------------------

event_list <- split(all_data, all_data$event_id)
n_total_events <- length(event_list)
cat("Fitting population models for", n_total_events, "events...\n")

# Parallel or sequential processing
if (opt$threads > 1 && .Platform$OS.type == "unix") {
    cat("Using", opt$threads, "threads\n")
    suppressPackageStartupMessages(library(parallel))
    outlier_list <- mclapply(
        event_list,
        fit_event_outliers,
        min_samples = opt$`min-samples`,
        zscore_threshold = opt$`zscore-threshold`,
        mc.cores = opt$threads
    )
} else {
    # Sequential with progress reporting
    outlier_list <- vector("list", n_total_events)
    names(outlier_list) <- names(event_list)
    progress_step <- max(1, n_total_events %/% 20)

    for (idx in seq_along(event_list)) {
        outlier_list[[idx]] <- fit_event_outliers(
            event_list[[idx]],
            min_samples = opt$`min-samples`,
            zscore_threshold = opt$`zscore-threshold`
        )
        if (idx %% progress_step == 0) {
            cat(sprintf("  Progress: %d/%d events (%.0f%%)\n",
                        idx, n_total_events, 100 * idx / n_total_events))
        }
    }
}

# Combine results
outlier_list <- outlier_list[!sapply(outlier_list, is.null)]

if (length(outlier_list) == 0) {
    cat("No outlier events detected.\n")
    write_empty_outputs(opt$`output-dir`)
    quit(save = "no", status = 0)
}

outlier_df <- do.call(rbind, outlier_list)
rownames(outlier_df) <- NULL

cat("Candidate outlier observations (pre-FDR):", nrow(outlier_df), "\n")

# ---- FDR correction ---------------------------------------------------------

outlier_df$padj <- p.adjust(outlier_df$pvalue, method = "BH")

# Filter by FDR threshold
sig_outliers <- outlier_df[outlier_df$padj < opt$`fdr-threshold`, ]

cat("Significant outliers (FDR <", opt$`fdr-threshold`, "):", nrow(sig_outliers), "\n")

# Sort by z-score magnitude
sig_outliers <- sig_outliers[order(-abs(sig_outliers$zscore)), ]

# Round p-values for readability
sig_outliers$pvalue <- signif(sig_outliers$pvalue, 4)
sig_outliers$padj <- signif(sig_outliers$padj, 4)

# ---- Write outlier_events.tsv -----------------------------------------------

output_events <- file.path(opt$`output-dir`, "outlier_events.tsv")
write.table(sig_outliers, output_events, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Outlier events written to:", output_events, "\n")

# ---- Generate per-sample outlier summary ------------------------------------

if (nrow(sig_outliers) > 0) {
    summary_df <- do.call(rbind, lapply(split(sig_outliers, sig_outliers$sample_id), function(sdf) {
        sample_id <- sdf$sample_id[1]
        n_outlier <- nrow(sdf)
        n_high <- sum(sdf$direction == "high_PSI")
        n_low <- sum(sdf$direction == "low_PSI")
        n_event_types <- length(unique(sdf$event_type))

        # Top outlier gene (by most extreme z-score)
        top_gene <- sdf$gene[which.max(abs(sdf$zscore))]

        data.frame(
            sample_id = sample_id,
            n_outlier_events = n_outlier,
            n_high_psi = n_high,
            n_low_psi = n_low,
            n_event_types_affected = n_event_types,
            top_outlier_gene = top_gene,
            stringsAsFactors = FALSE
        )
    }))

    # Sort by number of outlier events (descending)
    summary_df <- summary_df[order(-summary_df$n_outlier_events), ]
    rownames(summary_df) <- NULL
} else {
    summary_df <- data.frame(
        sample_id = character(),
        n_outlier_events = integer(),
        n_high_psi = integer(),
        n_low_psi = integer(),
        n_event_types_affected = integer(),
        top_outlier_gene = character(),
        stringsAsFactors = FALSE
    )
}

output_summary <- file.path(opt$`output-dir`, "outlier_summary.tsv")
write.table(summary_df, output_summary, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Sample outlier summary written to:", output_summary, "\n")

# ---- Print summary statistics -----------------------------------------------

cat("\n========================================\n")
cat("Bisbee Outlier Analysis Summary\n")
cat("========================================\n")
cat("Total events analyzed:", n_total_events, "\n")
cat("Events with outlier calls:", length(unique(sig_outliers$event_id)), "\n")
cat("Total outlier event-sample pairs:", nrow(sig_outliers), "\n")
cat("Samples with at least one outlier:", nrow(summary_df), "\n")

if (nrow(summary_df) > 0) {
    cat("\nTop 10 samples by outlier count:\n")
    top10 <- head(summary_df, 10)
    for (i in seq_len(nrow(top10))) {
        cat(sprintf("  %s: %d outliers (high: %d, low: %d) | top gene: %s\n",
                    top10$sample_id[i], top10$n_outlier_events[i],
                    top10$n_high_psi[i], top10$n_low_psi[i],
                    top10$top_outlier_gene[i]))
    }
}

# Event type breakdown
if (nrow(sig_outliers) > 0) {
    cat("\nOutlier counts by event type:\n")
    et_table <- table(sig_outliers$event_type)
    for (et in names(sort(et_table, decreasing = TRUE))) {
        cat(sprintf("  %s: %d\n", et, et_table[et]))
    }
}

cat("\nBisbee outlier detection complete.\n")
