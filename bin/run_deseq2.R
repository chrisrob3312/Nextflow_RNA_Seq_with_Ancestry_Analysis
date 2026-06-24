#!/usr/bin/env Rscript

# ============================================================================
# DESeq2 Differential Expression Analysis
# Supports: continuous ancestry, categorical contrasts, within-subgroup analyses
# Covariates: age, sex, blast_percentage, tumor_purity, batch, timepoint
# ============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(optparse)
    library(jsonlite)
    library(ggplot2)
    library(pheatmap)
    library(EnhancedVolcano)
    library(RColorBrewer)
    library(ggrepel)
})

# ---- Parse arguments ----
option_list <- list(
    make_option("--counts", type = "character", help = "Raw count matrix TSV"),
    make_option("--metadata", type = "character", help = "Sample metadata TSV"),
    make_option("--ancestry", type = "character", help = "Ancestry proportions TSV"),
    make_option("--contrasts", type = "character", help = "Contrasts JSON file"),
    make_option("--tumor-purity", type = "character", default = NULL,
                help = "TSV with sample_id and tumor_purity from ESTIMATE"),
    make_option("--covariates", type = "character", default = "batch,sex,age,tumor_purity",
                help = "Comma-separated covariate names for model"),
    make_option("--cytomolecular-subgroups", type = "character", default = "",
                help = "Comma-separated cytomolecular subgroups for within-group analysis"),
    make_option("--padj-threshold", type = "double", default = 0.05),
    make_option("--lfc-threshold", type = "double", default = 0.585),
    make_option("--output-dir", type = "character", default = "deseq2_results"),
    make_option("--plot-dir", type = "character", default = "deseq2_plots"),
    make_option("--rds-dir", type = "character", default = "deseq2_rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Load data ----
counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

# Load ancestry proportions and merge
if (!is.null(opt$ancestry) && file.exists(opt$ancestry)) {
    ancestry <- read.delim(opt$ancestry, check.names = FALSE)
    rownames(ancestry) <- ancestry$sample_id
    # Merge ancestry into metadata
    common_samples <- intersect(rownames(metadata), rownames(ancestry))
    ancestry_cols <- setdiff(colnames(ancestry), "sample_id")
    metadata[common_samples, ancestry_cols] <- ancestry[common_samples, ancestry_cols]
}

# Load ESTIMATE tumor purity and merge
if (!is.null(opt$`tumor-purity`) && file.exists(opt$`tumor-purity`)) {
    purity <- read.delim(opt$`tumor-purity`, check.names = FALSE)
    rownames(purity) <- purity$sample_id
    common_pur <- intersect(rownames(metadata), rownames(purity))
    metadata[common_pur, "tumor_purity"] <- purity[common_pur, "tumor_purity"]
    cat("Merged ESTIMATE tumor purity for", length(common_pur), "samples\n")
}

# Align samples
common_samples <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_samples]
metadata <- metadata[common_samples, ]

# ---- Load contrasts ----
if (!is.null(opt$contrasts) && file.exists(opt$contrasts)) {
    contrasts <- fromJSON(opt$contrasts)
} else {
    # Default contrasts for B-ALL with ancestry
    contrasts <- list(
        list(name = "relapse_vs_no_relapse",
             variable = "relapse_status", type = "categorical",
             reference = "no_relapse", target = "relapse"),
        list(name = "adi_q1_vs_rest",
             variable = "adi_quartile", type = "categorical",
             reference = "rest", target = "Q1"),
        list(name = "pct_african_continuous",
             variable = "pct_african", type = "continuous"),
        list(name = "pct_amerindigenous_continuous",
             variable = "pct_amerindigenous", type = "continuous"),
        list(name = "pct_east_asian_continuous",
             variable = "pct_east_asian", type = "continuous"),
        list(name = "pct_south_asian_continuous",
             variable = "pct_south_asian", type = "continuous"),
        list(name = "pct_european_continuous",
             variable = "pct_european", type = "continuous"),
        list(name = "graf_category",
             variable = "graf_category", type = "categorical_multi",
             reference = "EUR")
    )
}

# ---- Parse covariates ----
covariate_names <- trimws(unlist(strsplit(opt$covariates, ",")))
# Only keep covariates that exist in metadata and have >1 level
valid_covariates <- sapply(covariate_names, function(cov) {
    if (cov %in% colnames(metadata)) {
        vals <- metadata[[cov]]
        vals <- vals[!is.na(vals) & vals != "NA"]
        return(length(unique(vals)) > 1)
    }
    return(FALSE)
})
covariate_names <- covariate_names[valid_covariates]

cat("Using covariates:", paste(covariate_names, collapse = ", "), "\n")

# ---- Helper: run DESeq2 for a given contrast ----
run_deseq2_contrast <- function(counts, metadata, contrast, covariates, subset_var = NULL, subset_val = NULL) {
    meta <- metadata
    cts <- counts

    # Subset if within-group analysis
    if (!is.null(subset_var) && !is.null(subset_val)) {
        keep <- meta[[subset_var]] == subset_val
        meta <- meta[keep, , drop = FALSE]
        cts <- cts[, rownames(meta), drop = FALSE]
    }

    # Remove samples with NA for the variable of interest
    var_name <- contrast$variable
    if (!var_name %in% colnames(meta)) {
        warning(paste("Variable", var_name, "not found in metadata. Skipping."))
        return(NULL)
    }

    keep <- !is.na(meta[[var_name]]) & meta[[var_name]] != "NA"
    meta <- meta[keep, , drop = FALSE]
    cts <- cts[, rownames(meta), drop = FALSE]

    if (nrow(meta) < 4) {
        warning(paste("Too few samples for contrast", contrast$name))
        return(NULL)
    }

    # Build design formula
    covar_in_model <- covariates[covariates %in% colnames(meta)]
    # Remove covariates with only 1 level in this subset
    covar_in_model <- sapply(covar_in_model, function(cov) {
        vals <- meta[[cov]]
        vals <- vals[!is.na(vals) & vals != "NA"]
        length(unique(vals)) > 1
    })
    covar_in_model <- names(covar_in_model)[covar_in_model]

    if (contrast$type == "continuous") {
        meta[[var_name]] <- as.numeric(meta[[var_name]])
        formula_str <- paste("~", paste(c(covar_in_model, var_name), collapse = " + "))
    } else if (contrast$type == "categorical") {
        # For binary categorical: set up as factor
        if (!is.null(contrast$target) && !is.null(contrast$reference)) {
            if (contrast$reference == "rest") {
                meta[[var_name]] <- ifelse(meta[[var_name]] == contrast$target, contrast$target, "rest")
            }
            meta[[var_name]] <- factor(meta[[var_name]], levels = c(contrast$reference, contrast$target))
        } else {
            meta[[var_name]] <- as.factor(meta[[var_name]])
        }
        formula_str <- paste("~", paste(c(covar_in_model, var_name), collapse = " + "))
    } else if (contrast$type == "categorical_multi") {
        meta[[var_name]] <- as.factor(meta[[var_name]])
        if (!is.null(contrast$reference)) {
            meta[[var_name]] <- relevel(meta[[var_name]], ref = contrast$reference)
        }
        formula_str <- paste("~", paste(c(covar_in_model, var_name), collapse = " + "))
    }

    # Convert covariates to appropriate types
    numeric_covariates <- c("age", "tumor_purity", "blast_percentage")
    for (cov in covar_in_model) {
        if (cov %in% numeric_covariates) {
            meta[[cov]] <- as.numeric(meta[[cov]])
        } else {
            meta[[cov]] <- as.factor(meta[[cov]])
        }
    }

    # Filter genes with low counts
    keep_genes <- rowSums(cts >= 10) >= max(3, ncol(cts) * 0.1)
    cts <- cts[keep_genes, ]

    # Create DESeqDataSet
    design <- as.formula(formula_str)
    dds <- DESeqDataSetFromMatrix(countData = round(cts),
                                   colData = meta,
                                   design = design)

    # Run DESeq2
    dds <- DESeq(dds)

    # Extract results
    if (contrast$type == "continuous") {
        res <- results(dds, name = var_name, alpha = opt$`padj-threshold`)
    } else if (contrast$type == "categorical") {
        coef_name <- paste0(var_name, "_", contrast$target, "_vs_", contrast$reference)
        res <- results(dds, name = coef_name, alpha = opt$`padj-threshold`)
    } else if (contrast$type == "categorical_multi") {
        # Return results for all levels vs reference
        result_names <- resultsNames(dds)
        var_results <- grep(var_name, result_names, value = TRUE)
        res_list <- lapply(var_results, function(rn) {
            results(dds, name = rn, alpha = opt$`padj-threshold`)
        })
        names(res_list) <- var_results
        return(list(dds = dds, results = res_list, multi = TRUE))
    }

    return(list(dds = dds, results = res, multi = FALSE))
}

# ---- Run all contrasts ----
all_results <- list()
dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`rds-dir`, recursive = TRUE, showWarnings = FALSE)

# Summary tracking
summary_df <- data.frame(
    contrast = character(),
    subgroup = character(),
    n_samples = integer(),
    n_sig_up = integer(),
    n_sig_down = integer(),
    n_sig_total = integer(),
    stringsAsFactors = FALSE
)

for (i in seq_along(contrasts)) {
    contrast <- contrasts[[i]]
    cat("\n==== Running contrast:", contrast$name, "====\n")

    # Handle subset_variable / subset_value for stratified analyses (e.g., MRD-neg only)
    subset_var <- contrast$subset_variable
    subset_val <- contrast$subset_value
    result <- run_deseq2_contrast(counts, metadata, contrast, covariate_names,
                                  subset_var = subset_var, subset_val = subset_val)

    if (is.null(result)) next

    if (result$multi) {
        for (rn in names(result$results)) {
            res_df <- as.data.frame(result$results[[rn]])
            res_df <- res_df[order(res_df$padj), ]
            out_name <- paste0(contrast$name, "_", gsub("[^a-zA-Z0-9]", "_", rn))
            write.table(res_df, file.path(opt$`output-dir`, paste0(out_name, ".tsv")),
                        sep = "\t", quote = FALSE)

            n_up <- sum(res_df$padj < opt$`padj-threshold` & res_df$log2FoldChange > opt$`lfc-threshold`, na.rm = TRUE)
            n_down <- sum(res_df$padj < opt$`padj-threshold` & res_df$log2FoldChange < -opt$`lfc-threshold`, na.rm = TRUE)
            summary_df <- rbind(summary_df, data.frame(
                contrast = out_name, subgroup = "all",
                n_samples = ncol(counts), n_sig_up = n_up,
                n_sig_down = n_down, n_sig_total = n_up + n_down
            ))
        }
    } else {
        res_df <- as.data.frame(result$results)
        res_df <- res_df[order(res_df$padj), ]
        write.table(res_df, file.path(opt$`output-dir`, paste0(contrast$name, ".tsv")),
                    sep = "\t", quote = FALSE)

        # Save RDS
        saveRDS(result$dds, file.path(opt$`rds-dir`, paste0(contrast$name, "_dds.rds")))

        n_up <- sum(res_df$padj < opt$`padj-threshold` & res_df$log2FoldChange > opt$`lfc-threshold`, na.rm = TRUE)
        n_down <- sum(res_df$padj < opt$`padj-threshold` & res_df$log2FoldChange < -opt$`lfc-threshold`, na.rm = TRUE)
        summary_df <- rbind(summary_df, data.frame(
            contrast = contrast$name, subgroup = "all",
            n_samples = ncol(counts), n_sig_up = n_up,
            n_sig_down = n_down, n_sig_total = n_up + n_down
        ))

        # Volcano plot
        tryCatch({
            pdf(file.path(opt$`plot-dir`, paste0(contrast$name, "_volcano.pdf")), width = 10, height = 8)
            print(EnhancedVolcano(res_df,
                lab = rownames(res_df),
                x = 'log2FoldChange', y = 'padj',
                pCutoff = opt$`padj-threshold`,
                FCcutoff = opt$`lfc-threshold`,
                title = contrast$name,
                subtitle = paste0("Up: ", n_up, " | Down: ", n_down)))
            dev.off()
        }, error = function(e) { cat("Volcano plot failed:", e$message, "\n") })

        # MA plot
        tryCatch({
            pdf(file.path(opt$`plot-dir`, paste0(contrast$name, "_MA.pdf")), width = 10, height = 6)
            plotMA(result$results, main = contrast$name, ylim = c(-5, 5))
            dev.off()
        }, error = function(e) { cat("MA plot failed:", e$message, "\n") })
    }

    # ---- Within cytomolecular subgroup analyses ----
    if (opt$`cytomolecular-subgroups` != "") {
        subgroups <- trimws(unlist(strsplit(opt$`cytomolecular-subgroups`, ",")))
        for (sg in subgroups) {
            cat("  Running within subgroup:", sg, "\n")
            sg_result <- run_deseq2_contrast(counts, metadata, contrast, covariate_names,
                                              subset_var = "cytomolecular_subgroup", subset_val = sg)
            if (!is.null(sg_result) && !sg_result$multi) {
                sg_res_df <- as.data.frame(sg_result$results)
                sg_res_df <- sg_res_df[order(sg_res_df$padj), ]
                sg_name <- paste0(contrast$name, "_within_", gsub("[^a-zA-Z0-9]", "_", sg))
                write.table(sg_res_df, file.path(opt$`output-dir`, paste0(sg_name, ".tsv")),
                            sep = "\t", quote = FALSE)

                n_up <- sum(sg_res_df$padj < opt$`padj-threshold` & sg_res_df$log2FoldChange > opt$`lfc-threshold`, na.rm = TRUE)
                n_down <- sum(sg_res_df$padj < opt$`padj-threshold` & sg_res_df$log2FoldChange < -opt$`lfc-threshold`, na.rm = TRUE)
                summary_df <- rbind(summary_df, data.frame(
                    contrast = contrast$name, subgroup = sg,
                    n_samples = nrow(sg_result$dds@colData),
                    n_sig_up = n_up, n_sig_down = n_down, n_sig_total = n_up + n_down
                ))
            }
        }
    }
}

# Write summary
write.table(summary_df, file.path(opt$`output-dir`, "de_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDESeq2 analysis complete. Summary:\n")
print(summary_df)
