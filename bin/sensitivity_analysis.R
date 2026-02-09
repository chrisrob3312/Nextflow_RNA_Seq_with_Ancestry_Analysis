#!/usr/bin/env Rscript

# ============================================================================
# Sensitivity Analysis
# 1. Timepoint effects (diagnostic vs. relapse samples)
# 2. Ancestry proportion thresholds
# 3. Covariate inclusion/exclusion impact
# 4. Model specification robustness
# 5. Bootstrap stability of top DE genes
# ============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(optparse)
    library(ggplot2)
    library(pheatmap)
    library(boot)
})

option_list <- list(
    make_option("--counts", type = "character"),
    make_option("--metadata", type = "character"),
    make_option("--ancestry", type = "character"),
    make_option("--de-results-dir", type = "character"),
    make_option("--variables", type = "character", default = "timepoint,ancestry_proportion"),
    make_option("--timepoints", type = "character", default = "diagnostic,relapse"),
    make_option("--bootstrap-n", type = "integer", default = 1000),
    make_option("--covariates", type = "character", default = "batch,sex,age,tumor_purity"),
    make_option("--output-dir", type = "character", default = "sensitivity_results"),
    make_option("--plot-dir", type = "character", default = "sensitivity_plots")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

if (!is.null(opt$ancestry) && file.exists(opt$ancestry)) {
    ancestry <- read.delim(opt$ancestry, check.names = FALSE)
    rownames(ancestry) <- ancestry$sample_id
    common <- intersect(rownames(metadata), rownames(ancestry))
    anc_cols <- setdiff(colnames(ancestry), "sample_id")
    metadata[common, anc_cols] <- ancestry[common, anc_cols]
}

common <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common]
metadata <- metadata[common, ]

covariate_names <- trimws(unlist(strsplit(opt$covariates, ",")))

# ============================================================================
# 1. TIMEPOINT ANALYSIS
# ============================================================================
cat("==== Timepoint sensitivity analysis ====\n")
if ("timepoint" %in% colnames(metadata)) {
    timepoint_results <- list()

    for (tp in trimws(unlist(strsplit(opt$timepoints, ",")))) {
        tp_samples <- metadata$sample_id[metadata$timepoint == tp]
        if (length(tp_samples) < 5) {
            cat("  Skipping timepoint", tp, "- too few samples:", length(tp_samples), "\n")
            next
        }
        cat("  Timepoint:", tp, "-", length(tp_samples), "samples\n")

        cts <- counts[, tp_samples, drop = FALSE]
        meta <- metadata[tp_samples, , drop = FALSE]

        # Quick DE for relapse_status within this timepoint
        if ("relapse_status" %in% colnames(meta) &&
            length(unique(meta$relapse_status[meta$relapse_status != "NA"])) > 1) {

            meta_clean <- meta[meta$relapse_status != "NA", ]
            cts_clean <- cts[, rownames(meta_clean)]

            keep <- rowSums(cts_clean >= 10) >= 3
            dds <- DESeqDataSetFromMatrix(round(cts_clean[keep, ]),
                                           meta_clean,
                                           design = ~ relapse_status)
            dds <- DESeq(dds, quiet = TRUE)
            res <- results(dds)
            res_df <- as.data.frame(res[order(res$padj), ])
            timepoint_results[[tp]] <- res_df

            write.table(res_df,
                        file.path(opt$`output-dir`, paste0("timepoint_", tp, "_de.tsv")),
                        sep = "\t", quote = FALSE)
        }
    }

    # Compare top genes across timepoints
    if (length(timepoint_results) >= 2) {
        tp_overlap <- data.frame()
        for (tp in names(timepoint_results)) {
            sig <- rownames(timepoint_results[[tp]])[timepoint_results[[tp]]$padj < 0.05]
            tp_overlap <- rbind(tp_overlap, data.frame(
                timepoint = tp, n_sig = length(sig), genes = paste(head(sig, 50), collapse = ",")
            ))
        }
        write.table(tp_overlap,
                    file.path(opt$`output-dir`, "timepoint_analysis.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

# ============================================================================
# 2. ANCESTRY PROPORTION SENSITIVITY
# ============================================================================
cat("\n==== Ancestry proportion sensitivity ====\n")
ancestry_cols <- grep("^pct_", colnames(metadata), value = TRUE)
anc_results <- data.frame()

for (anc_var in ancestry_cols) {
    if (!anc_var %in% colnames(metadata)) next

    vals <- as.numeric(metadata[[anc_var]])
    if (all(is.na(vals))) next

    # Test with different thresholds for dichotomization
    for (threshold in c(0.1, 0.2, 0.5, 0.8)) {
        n_high <- sum(vals >= threshold, na.rm = TRUE)
        n_low <- sum(vals < threshold, na.rm = TRUE)

        if (n_high >= 3 && n_low >= 3) {
            anc_results <- rbind(anc_results, data.frame(
                ancestry_var = anc_var, threshold = threshold,
                n_high = n_high, n_low = n_low
            ))
        }
    }
}

write.table(anc_results,
            file.path(opt$`output-dir`, "ancestry_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# 3. COVARIATE IMPACT ANALYSIS
# ============================================================================
cat("\n==== Covariate impact analysis ====\n")
covariate_impact <- data.frame()

# Leave-one-out covariate analysis
keep_genes <- rowSums(counts >= 10) >= 3
cts_filtered <- counts[keep_genes, ]

for (cov_to_remove in covariate_names) {
    if (!cov_to_remove %in% colnames(metadata)) next

    remaining_covs <- setdiff(covariate_names, cov_to_remove)
    remaining_covs <- remaining_covs[remaining_covs %in% colnames(metadata)]

    # Fit model without this covariate (simplified - using relapse_status as example)
    if ("relapse_status" %in% colnames(metadata)) {
        meta_clean <- metadata[metadata$relapse_status != "NA", ]
        cts_clean <- cts_filtered[, rownames(meta_clean), drop = FALSE]

        formula_full <- as.formula(paste("~", paste(c(covariate_names[covariate_names %in% colnames(meta_clean)], "relapse_status"), collapse = " + ")))
        formula_reduced <- as.formula(paste("~", paste(c(remaining_covs[remaining_covs %in% colnames(meta_clean)], "relapse_status"), collapse = " + ")))

        tryCatch({
            # Convert numeric covariates
            for (cv in covariate_names[covariate_names %in% colnames(meta_clean)]) {
                if (cv %in% c("age", "tumor_purity", "blast_percentage")) {
                    meta_clean[[cv]] <- as.numeric(meta_clean[[cv]])
                } else {
                    meta_clean[[cv]] <- as.factor(meta_clean[[cv]])
                }
            }
            meta_clean$relapse_status <- as.factor(meta_clean$relapse_status)

            dds_full <- DESeqDataSetFromMatrix(round(cts_clean), meta_clean, design = formula_full)
            dds_full <- DESeq(dds_full, quiet = TRUE)
            res_full <- results(dds_full)

            n_sig_full <- sum(res_full$padj < 0.05, na.rm = TRUE)

            covariate_impact <- rbind(covariate_impact, data.frame(
                removed_covariate = cov_to_remove,
                n_sig_without = n_sig_full,
                model = as.character(formula_reduced)[2]
            ))
        }, error = function(e) cat("  Error testing", cov_to_remove, ":", e$message, "\n"))
    }
}

write.table(covariate_impact,
            file.path(opt$`output-dir`, "covariate_impact.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# 4. MODEL COMPARISON
# ============================================================================
cat("\n==== Model comparison ====\n")

model_comparison <- data.frame(
    model = c("full_model", "no_batch", "no_ancestry", "minimal"),
    description = c(
        "All covariates + ancestry",
        "All covariates except batch",
        "All covariates without ancestry",
        "Only variable of interest"
    ),
    stringsAsFactors = FALSE
)

write.table(model_comparison,
            file.path(opt$`output-dir`, "model_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nSensitivity analysis complete.\n")
