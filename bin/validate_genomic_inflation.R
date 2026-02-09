#!/usr/bin/env Rscript

# ============================================================================
# Genomic Inflation Factor Validation
# Checks for p-value inflation/deflation in differential expression results:
#   - Computes lambda (genomic inflation factor) from p-value distribution
#   - Generates QQ plot with 95% confidence band
#   - Flags lambda > 1.5 (potential confounding) or < 0.5 (overcorrection)
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
})

# ---- Parse arguments ----
option_list <- list(
    make_option("--de-results", type = "character",
                help = "DE results TSV file (must contain a 'pvalue' or 'P.Value' column)"),
    make_option("--contrast", type = "character", default = "unknown",
                help = "Contrast name for labeling outputs [default: %default]"),
    make_option("--output-dir", type = "character", default = "inflation_results",
                help = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate inputs ----
if (is.null(opt$`de-results`)) {
    stop("--de-results argument is required.")
}

if (!file.exists(opt$`de-results`)) {
    stop(paste("DE results file not found:", opt$`de-results`))
}

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- Load DE results ----
cat("Loading DE results:", opt$`de-results`, "\n")
de_res <- read.delim(opt$`de-results`, check.names = FALSE, row.names = 1)

# Find p-value column
pval_col <- NULL
for (candidate in c("pvalue", "P.Value", "PValue", "p.value", "p_value", "pval")) {
    if (candidate %in% colnames(de_res)) {
        pval_col <- candidate
        break
    }
}

if (is.null(pval_col)) {
    stop(paste("No p-value column found. Available columns:",
               paste(colnames(de_res), collapse = ", "),
               "\nExpected one of: pvalue, P.Value, PValue, p.value, p_value, pval"))
}

cat("Using p-value column:", pval_col, "\n")

# Extract and clean p-values
pvalues <- as.numeric(de_res[[pval_col]])
n_total <- length(pvalues)
n_na <- sum(is.na(pvalues))
pvalues <- pvalues[!is.na(pvalues)]

# Remove exact 0 and 1 p-values (can cause issues with qchisq)
n_zero <- sum(pvalues == 0)
n_one <- sum(pvalues == 1)
pvalues_clean <- pvalues[pvalues > 0 & pvalues < 1]

cat("Total genes:", n_total, "\n")
cat("  NA p-values:", n_na, "\n")
cat("  Zero p-values:", n_zero, "\n")
cat("  P-values used:", length(pvalues_clean), "\n")

if (length(pvalues_clean) < 10) {
    findings <- c(
        paste0("# Genomic Inflation Validation: ", opt$contrast),
        "",
        "## Status: SKIPPED",
        "",
        paste0("Too few valid p-values (n=", length(pvalues_clean),
               ") for genomic inflation analysis. At least 10 are required."),
        "",
        paste0("- Total genes: ", n_total),
        paste0("- NA p-values: ", n_na),
        paste0("- Zero p-values: ", n_zero),
        ""
    )
    writeLines(findings, file.path(opt$`output-dir`, "inflation_findings.md"))
    cat("Too few p-values. Findings written.\n")
    quit(status = 0)
}

# ---- Compute genomic inflation factor (lambda) ----
# Lambda is the ratio of the median observed chi-squared statistic to
# the expected median under the null (chi-squared with 1 df)
chisq_obs <- qchisq(1 - pvalues_clean, df = 1)
lambda <- median(chisq_obs) / qchisq(0.5, df = 1)

cat("Genomic inflation factor (lambda):", round(lambda, 4), "\n")

# Also compute lambda at the 10th percentile (lambda_10) for sensitivity
chisq_p10 <- quantile(chisq_obs, 0.9)
lambda_10 <- chisq_p10 / qchisq(0.1, df = 1)

# Mean chi-squared (alternative inflation measure)
mean_chisq <- mean(chisq_obs)

# ---- P-value distribution statistics ----
# Kolmogorov-Smirnov test for uniformity under null
ks_test <- tryCatch({
    ks.test(pvalues_clean, "punif")
}, error = function(e) {
    list(statistic = NA, p.value = NA)
})

# Proportion of p-values in bins
n_p <- length(pvalues_clean)
pval_bins <- c(
    p_lt_0.001 = sum(pvalues_clean < 0.001) / n_p * 100,
    p_lt_0.01 = sum(pvalues_clean < 0.01) / n_p * 100,
    p_lt_0.05 = sum(pvalues_clean < 0.05) / n_p * 100,
    p_gt_0.95 = sum(pvalues_clean > 0.95) / n_p * 100
)

# ---- Write inflation summary ----
summary_df <- data.frame(
    contrast = opt$contrast,
    n_genes_total = n_total,
    n_genes_tested = length(pvalues_clean),
    n_na_pvalues = n_na,
    n_zero_pvalues = n_zero,
    lambda = round(lambda, 4),
    lambda_10 = round(lambda_10, 4),
    mean_chisq = round(mean_chisq, 4),
    ks_statistic = round(as.numeric(ks_test$statistic), 4),
    ks_pvalue = signif(as.numeric(ks_test$p.value), 4),
    pct_p_lt_0.001 = round(pval_bins["p_lt_0.001"], 2),
    pct_p_lt_0.01 = round(pval_bins["p_lt_0.01"], 2),
    pct_p_lt_0.05 = round(pval_bins["p_lt_0.05"], 2),
    pct_p_gt_0.95 = round(pval_bins["p_gt_0.95"], 2),
    stringsAsFactors = FALSE
)

write.table(summary_df, file.path(opt$`output-dir`, "inflation_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary written to inflation_summary.tsv\n")

# ---- Generate QQ plot ----
tryCatch({
    pdf(file.path(opt$`output-dir`, "qq_plot.pdf"), width = 8, height = 8)

    n <- length(pvalues_clean)
    observed <- -log10(sort(pvalues_clean))
    expected <- -log10(ppoints(n))

    # 95% confidence band under the null
    # Based on order statistics of the uniform distribution
    upper_ci <- -log10(qbeta(0.025, 1:n, n:1))
    lower_ci <- -log10(qbeta(0.975, 1:n, n:1))

    # Determine axis limits
    max_val <- max(c(observed, expected, upper_ci), na.rm = TRUE)
    max_val <- min(max_val, 50)  # Cap at 50 for readability

    qq_df <- data.frame(
        expected = expected,
        observed = pmin(observed, 50),
        upper_ci = pmin(upper_ci, 50),
        lower_ci = pmin(lower_ci, 50)
    )

    p <- ggplot(qq_df, aes(x = expected, y = observed)) +
        geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
                     fill = "grey80", alpha = 0.5) +
        geom_point(size = 0.8, alpha = 0.5, color = "black") +
        geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
        labs(
            title = paste0("QQ Plot: ", opt$contrast),
            subtitle = paste0("lambda = ", round(lambda, 3),
                             " | n = ", n,
                             ifelse(lambda > 1.5, " | WARNING: inflated",
                                    ifelse(lambda < 0.5, " | WARNING: deflated", ""))),
            x = expression(Expected ~ -log[10](p)),
            y = expression(Observed ~ -log[10](p))
        ) +
        theme_bw(base_size = 12) +
        coord_cartesian(xlim = c(0, max(expected) * 1.05),
                        ylim = c(0, max_val * 1.05))

    print(p)

    # Also generate a p-value histogram
    hist_df <- data.frame(pvalue = pvalues_clean)
    p_hist <- ggplot(hist_df, aes(x = pvalue)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white", boundary = 0) +
        geom_hline(yintercept = n / 50, color = "red", linetype = "dashed") +
        labs(
            title = paste0("P-value Distribution: ", opt$contrast),
            subtitle = paste0("Red line = expected under null | lambda = ", round(lambda, 3)),
            x = "P-value",
            y = "Count"
        ) +
        theme_bw(base_size = 12)

    print(p_hist)

    dev.off()
    cat("QQ plot and p-value histogram written to qq_plot.pdf\n")
}, error = function(e) {
    cat("QQ plot generation failed:", e$message, "\n")
})

# ---- Generate findings markdown ----
findings <- c(
    paste0("# Genomic Inflation Validation: ", opt$contrast),
    "",
    "## Summary",
    "",
    paste0("| Metric | Value |"),
    paste0("|--------|-------|"),
    paste0("| Genes tested | ", length(pvalues_clean), " |"),
    paste0("| Lambda (genomic inflation factor) | ", round(lambda, 4), " |"),
    paste0("| Lambda_10 (10th percentile) | ", round(lambda_10, 4), " |"),
    paste0("| KS test statistic | ", round(as.numeric(ks_test$statistic), 4), " |"),
    paste0("| KS test p-value | ", signif(as.numeric(ks_test$p.value), 4), " |"),
    "",
    "## P-value Distribution",
    "",
    paste0("| P-value Range | Observed % | Expected under null % |"),
    paste0("|--------------|------------|----------------------|"),
    paste0("| p < 0.001 | ", round(pval_bins["p_lt_0.001"], 2), "% | 0.10% |"),
    paste0("| p < 0.01 | ", round(pval_bins["p_lt_0.01"], 2), "% | 1.00% |"),
    paste0("| p < 0.05 | ", round(pval_bins["p_lt_0.05"], 2), "% | 5.00% |"),
    paste0("| p > 0.95 | ", round(pval_bins["p_gt_0.95"], 2), "% | 5.00% |"),
    ""
)

# Warnings
warnings_list <- character(0)

if (lambda > 1.5) {
    warnings_list <- c(warnings_list,
        paste0("- **P-VALUE INFLATION** (lambda = ", round(lambda, 3),
               "): The genomic inflation factor exceeds 1.5, indicating systematic ",
               "p-value inflation. This may be caused by:"),
        "  - Unaccounted population structure or ancestry",
        "  - Batch effects not included in the model",
        "  - Sample relatedness",
        "  - Incorrect model specification",
        "",
        "  **Recommended actions:**",
        "  - Verify that ancestry covariates are included in the model",
        "  - Check for additional batch effects (sequencing center, library prep date)",
        "  - Consider using surrogate variable analysis (SVA) to capture hidden confounders",
        "  - Review the study design for unbalanced covariates",
        "")
} else if (lambda > 1.2) {
    warnings_list <- c(warnings_list,
        paste0("- **MILD INFLATION** (lambda = ", round(lambda, 3),
               "): Lambda is between 1.2 and 1.5. This represents mild inflation ",
               "that may reflect true biological signal or minor confounding. ",
               "Monitor but may not require correction."),
        "")
}

if (lambda < 0.5) {
    warnings_list <- c(warnings_list,
        paste0("- **P-VALUE DEFLATION** (lambda = ", round(lambda, 3),
               "): The genomic inflation factor is below 0.5, indicating possible ",
               "overcorrection. This may be caused by:"),
        "  - Too many covariates relative to sample size",
        "  - Overcorrection for ancestry or batch",
        "  - Covariates collinear with the variable of interest",
        "",
        "  **Recommended actions:**",
        "  - Review covariate selection; consider reducing the model",
        "  - Check for collinearity between covariates and the test variable",
        "  - Compare results with a simpler model",
        "")
} else if (lambda < 0.8) {
    warnings_list <- c(warnings_list,
        paste0("- **MILD DEFLATION** (lambda = ", round(lambda, 3),
               "): Lambda is between 0.5 and 0.8. The model may be slightly ",
               "overcorrected. Review covariates for potential collinearity ",
               "with the variable of interest."),
        "")
}

if (n_zero > n_total * 0.01) {
    warnings_list <- c(warnings_list,
        paste0("- **EXCESSIVE ZERO P-VALUES**: ", n_zero, " genes (",
               round(n_zero / n_total * 100, 1),
               "%) have p-value = 0. This is unusual and may indicate numerical ",
               "issues in the DE analysis. These genes were excluded from lambda ",
               "computation."),
        "")
}

if (n_na > n_total * 0.3) {
    warnings_list <- c(warnings_list,
        paste0("- **HIGH NA RATE**: ", round(n_na / n_total * 100, 1),
               "% of genes have NA p-values. This may indicate excessive ",
               "independent filtering or convergence failures in the DE model."),
        "")
}

# Check for anti-conservative p-value enrichment at low end
if (pval_bins["p_lt_0.05"] > 20) {
    warnings_list <- c(warnings_list,
        paste0("- **ENRICHED SMALL P-VALUES**: ", round(pval_bins["p_lt_0.05"], 1),
               "% of genes have p < 0.05 (expected ~5% under null). ",
               "Combined with the lambda value, this suggests ",
               ifelse(lambda > 1.2, "systematic confounding.", "strong true biological signal.")),
        "")
}

if (length(warnings_list) > 0) {
    findings <- c(findings, "## Warnings", "", warnings_list)
} else {
    findings <- c(findings,
        "## Status: PASS",
        "",
        paste0("Lambda = ", round(lambda, 3), " is within acceptable range (0.5 - 1.5). "),
        "The p-value distribution is consistent with a well-calibrated model.",
        "No evidence of systematic inflation or overcorrection.",
        "")
}

# Interpretation guide
findings <- c(findings,
    "## Interpretation Guide",
    "",
    "- **Lambda ~ 1.0**: P-values follow expected null distribution (well-calibrated)",
    "- **Lambda > 1.0**: More small p-values than expected (inflation)",
    "  - 1.0-1.2: Acceptable, may reflect true signal",
    "  - 1.2-1.5: Mild inflation, monitor for confounding",
    "  - >1.5: Significant inflation, likely confounding present",
    "- **Lambda < 1.0**: Fewer small p-values than expected (deflation/overcorrection)",
    "  - 0.8-1.0: Acceptable, possibly conservative",
    "  - 0.5-0.8: Mild overcorrection, review model",
    "  - <0.5: Severe overcorrection, simplify model",
    "",
    "The QQ plot shows observed vs expected -log10(p-values) with a 95% confidence",
    "band. Points above the diagonal indicate inflation; points below indicate deflation.",
    ""
)

writeLines(findings, file.path(opt$`output-dir`, "inflation_findings.md"))
cat("Findings written to inflation_findings.md\n")

cat("\nGenomic inflation validation complete.\n")
