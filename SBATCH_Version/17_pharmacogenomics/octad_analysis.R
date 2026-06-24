#!/usr/bin/env Rscript
# ============================================================================
# 17 - OCTAD Drug Reversal Analysis
# ============================================================================
# Uses the OCTAD Bioconductor package to identify drugs that reverse disease
# gene expression signatures by querying the LINCS L1000 connectivity map.
# For each DE contrast, creates a disease signature from top up/down-regulated
# genes and scores drugs using the sRGES (summarized Reversal Gene Expression
# Score) method.
#
# Required packages: octad, octad.db, ggplot2, tidyverse, optparse
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(octad)
    library(octad.db)
    library(tidyverse)
    library(ggplot2)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
    make_option("--de_results_dir", type = "character",
                help = "Directory containing DE results (TSV files with log2FC and padj)"),
    make_option("--output_dir", type = "character",
                help = "Output directory"),
    make_option("--n_top_genes", type = "integer", default = 250,
                help = "Number of top up/down genes for disease signature [default: %default]"),
    make_option("--n_drugs", type = "integer", default = 50,
                help = "Number of top drug candidates to report [default: %default]"),
    make_option("--threads", type = "integer", default = 8,
                help = "Number of threads [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$de_results_dir) || is.null(opt$output_dir)) {
    stop("Required: --de_results_dir, --output_dir")
}

cat("=== OCTAD Drug Reversal Analysis ===\n")
cat("DE results dir:", opt$de_results_dir, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Top genes per direction:", opt$n_top_genes, "\n")
cat("Top drugs to report:", opt$n_drugs, "\n")
cat("Threads:", opt$threads, "\n\n")

# ============================================================================
# Create output directories
# ============================================================================
drug_dir <- file.path(opt$output_dir, "drug_candidates")
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(drug_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Load DE result files
# ============================================================================
cat("Loading DE results...\n")

de_files <- list.files(opt$de_results_dir, pattern = "\\.tsv$",
                       full.names = TRUE, recursive = TRUE)
cat(sprintf("  Found %d DE result files\n", length(de_files)))

if (length(de_files) == 0) {
    stop("No DE result files found in: ", opt$de_results_dir)
}

# ============================================================================
# Process each contrast
# ============================================================================
all_drug_results <- list()

for (de_file in de_files) {
    contrast_name <- gsub("\\.tsv$", "", basename(de_file))
    cat(sprintf("\n--- Processing contrast: %s ---\n", contrast_name))

    tryCatch({
        # Load DE results
        de_df <- read.table(de_file, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)

        # Check required columns
        if (!all(c("padj", "log2FoldChange") %in% colnames(de_df))) {
            cat(sprintf("  Skipping %s: missing padj or log2FoldChange columns\n",
                        contrast_name))
            next
        }

        # Ensure gene names are available
        gene_col <- colnames(de_df)[1]
        de_df$gene <- de_df[[gene_col]]

        # Filter significant genes and sort
        de_sig <- de_df %>%
            filter(!is.na(padj), !is.na(log2FoldChange)) %>%
            arrange(padj)

        # Create disease signature: top up-regulated and top down-regulated genes
        up_genes <- de_sig %>%
            filter(log2FoldChange > 0) %>%
            arrange(padj, desc(log2FoldChange)) %>%
            slice_head(n = opt$n_top_genes) %>%
            pull(gene)

        down_genes <- de_sig %>%
            filter(log2FoldChange < 0) %>%
            arrange(padj, log2FoldChange) %>%
            slice_head(n = opt$n_top_genes) %>%
            pull(gene)

        cat(sprintf("  Disease signature: %d up-regulated, %d down-regulated genes\n",
                    length(up_genes), length(down_genes)))

        if (length(up_genes) < 10 || length(down_genes) < 10) {
            cat(sprintf("  WARNING: Too few genes for %s, skipping\n", contrast_name))
            next
        }

        # -----------------------------------------------------------------
        # Run sRGES computation via OCTAD
        # -----------------------------------------------------------------
        cat("  Computing sRGES drug reversal scores...\n")

        # Compute the reversal scores using LINCS L1000
        sRGES_results <- runsRGES(
            dz_signature = de_sig %>%
                select(gene, log2FoldChange, padj) %>%
                rename(Symbol = gene, log2FC = log2FoldChange),
            output_path = file.path(drug_dir, paste0(contrast_name, "_sRGES_raw")),
            permutations = 10000
        )

        if (!is.null(sRGES_results) && nrow(sRGES_results) > 0) {
            # Rank by sRGES (most negative = strongest reversal)
            sRGES_ranked <- sRGES_results %>%
                arrange(sRGES) %>%
                mutate(rank = row_number())

            # Save top N drug candidates
            top_candidates <- head(sRGES_ranked, opt$n_drugs)
            output_file <- file.path(drug_dir,
                                     paste0(contrast_name, "_top_drug_candidates.tsv"))
            write.table(top_candidates, output_file,
                        sep = "\t", quote = FALSE, row.names = FALSE)

            # Save full results
            full_output <- file.path(drug_dir,
                                     paste0(contrast_name, "_all_sRGES_scores.tsv"))
            write.table(sRGES_ranked, full_output,
                        sep = "\t", quote = FALSE, row.names = FALSE)

            cat(sprintf("  Saved %d drug candidates (top %d highlighted)\n",
                        nrow(sRGES_ranked), min(opt$n_drugs, nrow(sRGES_ranked))))

            # Store for combined output
            top_candidates$contrast <- contrast_name
            all_drug_results[[contrast_name]] <- top_candidates

            # -----------------------------------------------------------------
            # Visualization: bar plot of top candidates
            # -----------------------------------------------------------------
            cat("  Generating visualization...\n")

            # Determine drug name column
            drug_name_col <- intersect(c("pert_iname", "drug_name", "name"),
                                       colnames(top_candidates))
            if (length(drug_name_col) > 0) {
                drug_name_col <- drug_name_col[1]
            } else {
                drug_name_col <- colnames(top_candidates)[1]
            }

            plot_data <- head(top_candidates, 30)

            p <- ggplot(plot_data,
                        aes(x = reorder(!!sym(drug_name_col), -sRGES),
                            y = sRGES)) +
                geom_bar(stat = "identity",
                         fill = ifelse(plot_data$sRGES < 0, "steelblue", "coral")) +
                coord_flip() +
                theme_minimal() +
                theme(axis.text.y = element_text(size = 8)) +
                labs(
                    title = paste("Top Drug Candidates -", contrast_name),
                    subtitle = "Drugs ranked by sRGES (more negative = stronger reversal)",
                    x = "Drug",
                    y = "sRGES (Summarized Reversal Gene Expression Score)"
                ) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")

            ggsave(file.path(plot_dir,
                             paste0(contrast_name, "_top_drugs.pdf")),
                   p, width = 10, height = 8)

        } else {
            cat(sprintf("  WARNING: No sRGES results returned for %s\n",
                        contrast_name))
        }

    }, error = function(e) {
        cat(sprintf("  ERROR processing %s: %s\n", contrast_name, conditionMessage(e)))
    })
}

# ============================================================================
# Combined summary across contrasts
# ============================================================================
cat("\n--- Combined Summary ---\n")

if (length(all_drug_results) > 0) {
    combined <- bind_rows(all_drug_results)
    write.table(combined,
                file.path(drug_dir, "all_contrasts_top_candidates.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Combined results: %d drug candidates across %d contrasts\n",
                nrow(combined), length(all_drug_results)))

    # Summary: drugs appearing in multiple contrasts
    drug_name_col <- intersect(c("pert_iname", "drug_name", "name"),
                               colnames(combined))
    if (length(drug_name_col) > 0) {
        recurrent_drugs <- combined %>%
            count(!!sym(drug_name_col[1]), sort = TRUE) %>%
            filter(n > 1)

        if (nrow(recurrent_drugs) > 0) {
            write.table(recurrent_drugs,
                        file.path(drug_dir, "recurrent_drug_candidates.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            cat(sprintf("  %d drugs appear in multiple contrasts\n",
                        nrow(recurrent_drugs)))
        }
    }
} else {
    cat("  WARNING: No drug results generated from any contrast\n")
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== OCTAD Analysis Complete ===\n")
cat(sprintf("Output directory: %s\n", opt$output_dir))
cat("Key outputs:\n")
cat(sprintf("  - Drug candidates: %s\n", drug_dir))
cat(sprintf("  - Plots: %s\n", plot_dir))
