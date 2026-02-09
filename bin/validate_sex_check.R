#!/usr/bin/env Rscript

# ============================================================================
# Sex Validation Check
# Validates reported sex in metadata against expression of sex-linked genes:
#   - XIST: high expression in females (X-inactivation)
#   - DDX3Y, USP9Y, UTY, KDM5D: Y-chromosome genes, expressed in males
# Flags mismatches that may indicate sample swaps or metadata errors.
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
})

# ---- Parse arguments ----
option_list <- list(
    make_option("--counts", type = "character",
                help = "Normalized or raw count matrix TSV (genes x samples)"),
    make_option("--metadata", type = "character",
                help = "Sample metadata TSV with 'sample_id' and 'sex' columns"),
    make_option("--output-dir", type = "character", default = "sex_check_results",
                help = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate inputs ----
if (is.null(opt$counts)) {
    stop("--counts argument is required.")
}

if (is.null(opt$metadata)) {
    stop("--metadata argument is required.")
}

if (!file.exists(opt$counts)) {
    stop(paste("Counts file not found:", opt$counts))
}

if (!file.exists(opt$metadata)) {
    stop(paste("Metadata file not found:", opt$metadata))
}

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
cat("Loading count matrix:", opt$counts, "\n")
counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)

cat("Loading metadata:", opt$metadata, "\n")
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

# Align samples
common_samples <- intersect(colnames(counts), rownames(metadata))
if (length(common_samples) == 0) {
    stop("No common samples between count matrix and metadata.")
}
counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

cat("Samples in analysis:", length(common_samples), "\n")

# ---- Check for sex column ----
sex_col <- NULL
for (candidate in c("sex", "Sex", "gender", "Gender", "SEX")) {
    if (candidate %in% colnames(metadata)) {
        sex_col <- candidate
        break
    }
}

if (is.null(sex_col)) {
    # Write findings noting the absence of sex information
    findings <- c(
        "# Sex Check Validation",
        "",
        "## Status: SKIPPED",
        "",
        "No sex column found in metadata (checked: sex, Sex, gender, Gender, SEX).",
        "Cannot validate reported sex against expression. Consider adding sex",
        "information to the metadata for quality control purposes.",
        ""
    )
    writeLines(findings, file.path(opt$`output-dir`, "sex_check_findings.md"))
    cat("No sex column in metadata. Findings written.\n")
    quit(status = 0)
}

cat("Using sex column:", sex_col, "\n")
reported_sex <- metadata[[sex_col]]
cat("Reported sex distribution:\n")
print(table(reported_sex, useNA = "ifany"))

# ---- Define sex-linked genes ----
# XIST: X-inactivation specific transcript (high in females)
xist_genes <- c("XIST")

# Y-chromosome genes reliably expressed in males
y_genes <- c("DDX3Y", "USP9Y", "UTY", "KDM5D")

# Also check common aliases / ENSEMBL IDs that might appear
xist_aliases <- c("XIST", "ENSG00000229807")
y_aliases <- c("DDX3Y", "USP9Y", "UTY", "KDM5D",
               "ENSG00000067048", "ENSG00000114374",
               "ENSG00000183878", "ENSG00000012817")

# ---- Find available sex-linked genes ----
find_genes <- function(gene_list, count_rownames) {
    found <- character(0)
    for (g in gene_list) {
        # Exact match
        if (g %in% count_rownames) {
            found <- c(found, g)
            next
        }
        # Check if gene name appears as part of rowname (e.g., ENSEMBL|SYMBOL format)
        partial <- grep(paste0("(^|\\|)", g, "($|\\|)"), count_rownames, value = TRUE)
        if (length(partial) > 0) {
            found <- c(found, partial[1])
        }
    }
    return(unique(found))
}

available_xist <- find_genes(xist_aliases, rownames(counts))
available_y <- find_genes(y_aliases, rownames(counts))

cat("Available XIST gene(s):", paste(available_xist, collapse = ", "),
    ifelse(length(available_xist) == 0, "(none found)", ""), "\n")
cat("Available Y-chromosome gene(s):", paste(available_y, collapse = ", "),
    ifelse(length(available_y) == 0, "(none found)", ""), "\n")

if (length(available_xist) == 0 && length(available_y) == 0) {
    findings <- c(
        "# Sex Check Validation",
        "",
        "## Status: SKIPPED",
        "",
        "No sex-linked genes found in count matrix.",
        "",
        "Searched for:",
        paste0("- XIST genes: ", paste(xist_aliases, collapse = ", ")),
        paste0("- Y-chromosome genes: ", paste(y_aliases, collapse = ", ")),
        "",
        "This may indicate the count matrix uses a different gene ID format.",
        "Consider re-running with a count matrix that includes gene symbols.",
        ""
    )
    writeLines(findings, file.path(opt$`output-dir`, "sex_check_findings.md"))
    cat("No sex-linked genes found. Findings written.\n")
    quit(status = 0)
}

# ---- Compute sex-specific expression scores ----
# Log2-transform counts (adding pseudocount)
log2_counts <- log2(as.matrix(counts) + 1)

# XIST score: mean log2 expression of XIST genes
if (length(available_xist) > 0) {
    xist_expr <- colMeans(log2_counts[available_xist, , drop = FALSE])
} else {
    xist_expr <- rep(0, ncol(counts))
    names(xist_expr) <- colnames(counts)
}

# Y-chromosome score: mean log2 expression of Y genes
if (length(available_y) > 0) {
    y_expr <- colMeans(log2_counts[available_y, , drop = FALSE])
} else {
    y_expr <- rep(0, ncol(counts))
    names(y_expr) <- colnames(counts)
}

# ---- Infer sex from expression ----
# Standardize scores for classification
xist_z <- if (sd(xist_expr) > 0) scale(xist_expr)[, 1] else rep(0, length(xist_expr))
y_z <- if (sd(y_expr) > 0) scale(y_expr)[, 1] else rep(0, length(y_expr))

# Classification logic:
#   Female: high XIST (z > 0) and low Y (z < 0), or high XIST if Y not available
#   Male: low XIST (z < 0) and high Y (z > 0), or high Y if XIST not available
#   Ambiguous: mixed signals
infer_sex <- function(xist_z_val, y_z_val, has_xist, has_y) {
    if (has_xist && has_y) {
        if (xist_z_val > 0 && y_z_val < 0) return("F")
        if (xist_z_val < 0 && y_z_val > 0) return("M")
        return("ambiguous")
    } else if (has_xist) {
        if (xist_z_val > 0.5) return("F")
        if (xist_z_val < -0.5) return("M")
        return("ambiguous")
    } else if (has_y) {
        if (y_z_val > 0.5) return("M")
        if (y_z_val < -0.5) return("F")
        return("ambiguous")
    }
    return("unknown")
}

has_xist <- length(available_xist) > 0
has_y <- length(available_y) > 0

inferred_sex <- sapply(seq_along(common_samples), function(i) {
    infer_sex(xist_z[i], y_z[i], has_xist, has_y)
})

# ---- Standardize reported sex for comparison ----
standardize_sex <- function(sex_values) {
    sex_upper <- toupper(trimws(as.character(sex_values)))
    result <- rep(NA_character_, length(sex_upper))
    result[sex_upper %in% c("F", "FEMALE", "XX")] <- "F"
    result[sex_upper %in% c("M", "MALE", "XY")] <- "M"
    return(result)
}

reported_std <- standardize_sex(reported_sex)

# ---- Build results table ----
results_df <- data.frame(
    sample_id = common_samples,
    reported_sex = as.character(reported_sex),
    reported_sex_std = reported_std,
    inferred_sex = inferred_sex,
    xist_expression = round(xist_expr, 4),
    y_chromosome_expression = round(y_expr, 4),
    xist_zscore = round(xist_z, 4),
    y_zscore = round(y_z, 4),
    stringsAsFactors = FALSE
)

# Flag mismatches
results_df$mismatch <- FALSE
results_df$mismatch[!is.na(results_df$reported_sex_std) &
                     results_df$inferred_sex != "ambiguous" &
                     results_df$inferred_sex != "unknown" &
                     results_df$reported_sex_std != results_df$inferred_sex] <- TRUE

n_mismatches <- sum(results_df$mismatch)
n_ambiguous <- sum(results_df$inferred_sex == "ambiguous")
n_checked <- sum(!is.na(results_df$reported_sex_std) &
                  results_df$inferred_sex != "ambiguous" &
                  results_df$inferred_sex != "unknown")

cat("\nSex check results:\n")
cat("  Samples checked:", n_checked, "\n")
cat("  Mismatches:", n_mismatches, "\n")
cat("  Ambiguous:", n_ambiguous, "\n")

write.table(results_df, file.path(opt$`output-dir`, "sex_check_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Results written to sex_check_results.tsv\n")

# ---- Generate plot ----
tryCatch({
    pdf(file.path(opt$`output-dir`, "sex_check_plot.pdf"), width = 8, height = 7)

    plot_df <- results_df
    plot_df$display_sex <- ifelse(is.na(plot_df$reported_sex_std),
                                   paste0(plot_df$reported_sex, " (unknown)"),
                                   plot_df$reported_sex_std)
    plot_df$status <- ifelse(plot_df$mismatch, "MISMATCH",
                              ifelse(plot_df$inferred_sex == "ambiguous", "Ambiguous", "Match"))

    p <- ggplot(plot_df, aes(x = xist_expression, y = y_chromosome_expression)) +
        geom_point(aes(color = display_sex, shape = status), size = 3, alpha = 0.7) +
        scale_shape_manual(values = c("Match" = 16, "MISMATCH" = 4, "Ambiguous" = 1)) +
        labs(
            title = "Sex Check: XIST vs Y-chromosome Expression",
            subtitle = paste0("Mismatches: ", n_mismatches,
                             " | Ambiguous: ", n_ambiguous,
                             " | Checked: ", n_checked),
            x = paste0("XIST Expression (log2 + 1)"),
            y = paste0("Y-chromosome Expression (log2 + 1, mean of ",
                       paste(available_y, collapse = ", "), ")"),
            color = "Reported Sex",
            shape = "Status"
        ) +
        theme_bw(base_size = 12)

    # Label mismatched samples
    if (n_mismatches > 0) {
        mismatch_df <- plot_df[plot_df$mismatch, ]
        p <- p + geom_text(data = mismatch_df,
                            aes(label = sample_id),
                            vjust = -1, size = 2.5, color = "red")
    }

    print(p)
    dev.off()
    cat("Plot written to sex_check_plot.pdf\n")
}, error = function(e) {
    cat("Plot generation failed:", e$message, "\n")
})

# ---- Generate findings markdown ----
findings <- c(
    "# Sex Check Validation",
    "",
    "## Summary",
    "",
    paste0("| Metric | Value |"),
    paste0("|--------|-------|"),
    paste0("| Total samples | ", length(common_samples), " |"),
    paste0("| Samples with reported sex | ", sum(!is.na(reported_std)), " |"),
    paste0("| Samples checked (non-ambiguous) | ", n_checked, " |"),
    paste0("| Mismatches detected | ", n_mismatches, " |"),
    paste0("| Ambiguous classification | ", n_ambiguous, " |"),
    paste0("| XIST genes available | ", paste(available_xist, collapse = ", "),
           ifelse(length(available_xist) == 0, "none", ""), " |"),
    paste0("| Y-chromosome genes available | ", paste(available_y, collapse = ", "),
           ifelse(length(available_y) == 0, "none", ""), " |"),
    ""
)

# Warnings
warnings_list <- character(0)

if (n_mismatches > 0) {
    mismatch_samples <- results_df[results_df$mismatch, ]
    warnings_list <- c(warnings_list,
        paste0("- **SEX MISMATCH DETECTED**: ", n_mismatches, " sample(s) have discordant ",
               "reported and expression-inferred sex. These may indicate sample swaps, ",
               "metadata errors, or unusual biology (e.g., sex chromosome aneuploidies)."),
        "")

    # Detail the mismatched samples
    warnings_list <- c(warnings_list, "  Mismatched samples:", "")
    for (i in seq_len(nrow(mismatch_samples))) {
        s <- mismatch_samples[i, ]
        warnings_list <- c(warnings_list,
            paste0("  - **", s$sample_id, "**: reported=", s$reported_sex_std,
                   ", inferred=", s$inferred_sex,
                   " (XIST=", round(s$xist_expression, 2),
                   ", Y=", round(s$y_chromosome_expression, 2), ")"))
    }
    warnings_list <- c(warnings_list, "")
}

if (n_ambiguous > sum(is.na(reported_std))) {
    warnings_list <- c(warnings_list,
        paste0("- **AMBIGUOUS SAMPLES**: ", n_ambiguous, " sample(s) could not be ",
               "confidently classified as male or female based on expression. ",
               "This may indicate partial sex chromosome loss or mixed samples."),
        "")
}

mismatch_rate <- if (n_checked > 0) n_mismatches / n_checked * 100 else 0
if (mismatch_rate > 5) {
    warnings_list <- c(warnings_list,
        paste0("- **HIGH MISMATCH RATE**: ", round(mismatch_rate, 1),
               "% of checked samples show sex mismatches. This exceeds the 5% threshold ",
               "and suggests systematic metadata issues."),
        "")
}

if (length(warnings_list) > 0) {
    findings <- c(findings, "## Warnings", "", warnings_list)
} else {
    findings <- c(findings,
        "## Status: PASS",
        "",
        "All samples with reported sex show consistent expression of sex-linked genes.",
        "No sample swaps or metadata errors detected based on sex-linked gene expression.",
        "")
}

# Methods
findings <- c(findings,
    "## Methods",
    "",
    "Sex was inferred from RNA-seq expression of sex-linked genes:",
    "",
    "- **XIST**: X-inactivation specific transcript, highly expressed in females",
    "- **Y-chromosome genes** (DDX3Y, USP9Y, UTY, KDM5D): expressed in males",
    "",
    "Expression values were log2-transformed and z-scored. Samples were classified as:",
    "",
    "- **Female**: high XIST (z > 0) and low Y-chromosome expression (z < 0)",
    "- **Male**: low XIST (z < 0) and high Y-chromosome expression (z > 0)",
    "- **Ambiguous**: mixed or intermediate expression patterns",
    "",
    "Mismatches are flagged when inferred sex differs from reported sex in metadata.",
    ""
)

writeLines(findings, file.path(opt$`output-dir`, "sex_check_findings.md"))
cat("Findings written to sex_check_findings.md\n")

cat("\nSex check validation complete.\n")
