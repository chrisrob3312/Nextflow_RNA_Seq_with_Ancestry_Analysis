#!/usr/bin/env Rscript

# ============================================================================
# limma-voom Differential Expression Analysis
# Parallel implementation to DESeq2 for concordance checking
# ============================================================================

suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
    library(optparse)
    library(jsonlite)
    library(ggplot2)
    library(EnhancedVolcano)
})

option_list <- list(
    make_option("--counts", type = "character", help = "Raw count matrix TSV"),
    make_option("--metadata", type = "character", help = "Sample metadata TSV"),
    make_option("--ancestry", type = "character", help = "Ancestry proportions TSV"),
    make_option("--contrasts", type = "character", help = "Contrasts JSON file"),
    make_option("--covariates", type = "character", default = "batch,sex,age,blast_percentage"),
    make_option("--cytomolecular-subgroups", type = "character", default = ""),
    make_option("--padj-threshold", type = "double", default = 0.05),
    make_option("--lfc-threshold", type = "double", default = 0.585),
    make_option("--output-dir", type = "character", default = "limma_results"),
    make_option("--plot-dir", type = "character", default = "limma_plots"),
    make_option("--rds-dir", type = "character", default = "limma_rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Load data ----
counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

if (!is.null(opt$ancestry) && file.exists(opt$ancestry)) {
    ancestry <- read.delim(opt$ancestry, check.names = FALSE)
    rownames(ancestry) <- ancestry$sample_id
    common <- intersect(rownames(metadata), rownames(ancestry))
    ancestry_cols <- setdiff(colnames(ancestry), "sample_id")
    metadata[common, ancestry_cols] <- ancestry[common, ancestry_cols]
}

common_samples <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_samples]
metadata <- metadata[common_samples, ]

# ---- Load contrasts ----
if (!is.null(opt$contrasts) && file.exists(opt$contrasts)) {
    contrasts_list <- fromJSON(opt$contrasts)
} else {
    contrasts_list <- list(
        list(name = "relapse_vs_no_relapse", variable = "relapse_status",
             type = "categorical", reference = "no_relapse", target = "relapse"),
        list(name = "pct_african_continuous", variable = "pct_african", type = "continuous"),
        list(name = "pct_amerindigenous_continuous", variable = "pct_amerindigenous", type = "continuous"),
        list(name = "adi_q1_vs_rest", variable = "adi_quartile",
             type = "categorical", reference = "rest", target = "Q1")
    )
}

covariate_names <- trimws(unlist(strsplit(opt$covariates, ",")))

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`rds-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- Run limma-voom per contrast ----
run_limma_contrast <- function(counts, metadata, contrast, covariates, subset_var = NULL, subset_val = NULL) {
    meta <- metadata
    cts <- counts

    if (!is.null(subset_var) && !is.null(subset_val)) {
        keep <- meta[[subset_var]] == subset_val
        meta <- meta[keep, , drop = FALSE]
        cts <- cts[, rownames(meta), drop = FALSE]
    }

    var_name <- contrast$variable
    if (!var_name %in% colnames(meta)) return(NULL)

    keep <- !is.na(meta[[var_name]]) & meta[[var_name]] != "NA"
    meta <- meta[keep, , drop = FALSE]
    cts <- cts[, rownames(meta), drop = FALSE]

    if (nrow(meta) < 4) return(NULL)

    # Create DGEList and filter
    dge <- DGEList(counts = round(cts))
    keep_genes <- filterByExpr(dge, min.count = 10, min.total.count = 15)
    dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")

    # Build design
    covar_in_model <- covariates[covariates %in% colnames(meta)]
    covar_in_model <- covar_in_model[sapply(covar_in_model, function(cov) {
        vals <- meta[[cov]][!is.na(meta[[cov]]) & meta[[cov]] != "NA"]
        length(unique(vals)) > 1
    })]

    for (cov in covar_in_model) {
        if (cov %in% c("age", "tumor_purity", "blast_percentage")) {
            meta[[cov]] <- as.numeric(meta[[cov]])
        } else {
            meta[[cov]] <- as.factor(meta[[cov]])
        }
    }

    if (contrast$type == "continuous") {
        meta[[var_name]] <- as.numeric(meta[[var_name]])
    } else if (contrast$type == "categorical") {
        if (!is.null(contrast$reference) && contrast$reference == "rest") {
            meta[[var_name]] <- ifelse(meta[[var_name]] == contrast$target, contrast$target, "rest")
        }
        meta[[var_name]] <- factor(meta[[var_name]], levels = c(contrast$reference, contrast$target))
    } else {
        meta[[var_name]] <- as.factor(meta[[var_name]])
        if (!is.null(contrast$reference)) {
            meta[[var_name]] <- relevel(meta[[var_name]], ref = contrast$reference)
        }
    }

    formula_str <- paste("~0 +", paste(c(covar_in_model, var_name), collapse = " + "))
    design <- model.matrix(as.formula(formula_str), data = meta)

    # voom transformation
    v <- voom(dge, design, plot = FALSE)

    # Fit
    fit <- lmFit(v, design)
    fit <- eBayes(fit)

    return(list(fit = fit, voom = v, design = design))
}

summary_df <- data.frame(contrast = character(), subgroup = character(),
                          n_samples = integer(), n_sig_up = integer(),
                          n_sig_down = integer(), n_sig_total = integer(),
                          stringsAsFactors = FALSE)

for (i in seq_along(contrasts_list)) {
    contrast <- contrasts_list[[i]]
    cat("\n==== Running limma contrast:", contrast$name, "====\n")

    subset_var <- contrast$subset_variable
    subset_val <- contrast$subset_value
    result <- run_limma_contrast(counts, metadata, contrast, covariate_names,
                                 subset_var = subset_var, subset_val = subset_val)
    if (is.null(result)) next

    tt <- topTable(result$fit, number = Inf, sort.by = "P")
    write.table(tt, file.path(opt$`output-dir`, paste0(contrast$name, "_limma.tsv")),
                sep = "\t", quote = FALSE)
    saveRDS(result$fit, file.path(opt$`rds-dir`, paste0(contrast$name, "_fit.rds")))

    n_up <- sum(tt$adj.P.Val < opt$`padj-threshold` & tt$logFC > opt$`lfc-threshold`, na.rm = TRUE)
    n_down <- sum(tt$adj.P.Val < opt$`padj-threshold` & tt$logFC < -opt$`lfc-threshold`, na.rm = TRUE)
    summary_df <- rbind(summary_df, data.frame(
        contrast = contrast$name, subgroup = "all",
        n_samples = ncol(counts), n_sig_up = n_up,
        n_sig_down = n_down, n_sig_total = n_up + n_down
    ))

    # Volcano plot
    tryCatch({
        pdf(file.path(opt$`plot-dir`, paste0(contrast$name, "_limma_volcano.pdf")), width = 10, height = 8)
        print(EnhancedVolcano(tt, lab = rownames(tt),
            x = 'logFC', y = 'adj.P.Val',
            pCutoff = opt$`padj-threshold`, FCcutoff = opt$`lfc-threshold`,
            title = paste(contrast$name, "(limma-voom)")))
        dev.off()
    }, error = function(e) cat("Volcano plot error:", e$message, "\n"))
}

write.table(summary_df, file.path(opt$`output-dir`, "limma_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nlimma-voom analysis complete.\n")
