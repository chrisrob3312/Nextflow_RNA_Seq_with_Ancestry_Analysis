#!/usr/bin/env Rscript
# ============================================================================
# 17 - signatureSearch Drug Connectivity Analysis
# ============================================================================
# Uses signatureSearch and signatureSearchData Bioconductor packages to:
#   1. For each DE contrast, extract top 150 up + 150 down genes as query
#   2. Search the LINCS L1000 database using the CMAP connectivity method
#   3. Perform Drug Set Enrichment Analysis (DSEA) by mechanism of action
#   4. Save ranked drug candidates and MOA enrichment results
#
# Required packages: signatureSearch, signatureSearchData, ExperimentHub,
#                    ggplot2, tidyverse, optparse
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(signatureSearch)
    library(signatureSearchData)
    library(ExperimentHub)
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
    make_option("--n_up", type = "integer", default = 150,
                help = "Number of top up-regulated genes for query [default: %default]"),
    make_option("--n_down", type = "integer", default = 150,
                help = "Number of top down-regulated genes for query [default: %default]"),
    make_option("--threads", type = "integer", default = 8,
                help = "Number of threads [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$de_results_dir) || is.null(opt$output_dir)) {
    stop("Required: --de_results_dir, --output_dir")
}

cat("=== signatureSearch Drug Connectivity Analysis ===\n")
cat("DE results dir:", opt$de_results_dir, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Query size: top", opt$n_up, "up +", opt$n_down, "down genes\n")
cat("Threads:", opt$threads, "\n\n")

# ============================================================================
# Create output directories
# ============================================================================
ranking_dir <- file.path(opt$output_dir, "drug_rankings")
dsea_dir <- file.path(opt$output_dir, "dsea")
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(ranking_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dsea_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Load LINCS L1000 reference database
# ============================================================================
cat("Loading LINCS L1000 reference database...\n")

# Access the LINCS L1000 database via ExperimentHub
eh <- ExperimentHub()

# Load the LINCS L1000 signature database (level 5 - z-scores)
# signatureSearchData provides pre-built HDF5-backed databases
tryCatch({
    db_path <- system.file("extdata", "sample_db.h5", package = "signatureSearchData")
    if (db_path == "" || !file.exists(db_path)) {
        cat("  Downloading LINCS L1000 database from ExperimentHub...\n")
        lincs_db <- eh[["EH3228"]]  # LINCS L1000 Level 5 (landmark genes)
    } else {
        lincs_db <- db_path
    }
    cat("  LINCS L1000 database loaded.\n")
}, error = function(e) {
    cat(sprintf("  WARNING: Could not load pre-built DB, building from ExperimentHub: %s\n",
                conditionMessage(e)))
    lincs_db <- NULL
})

# Build the SignatureSearch database object
cat("  Building search database...\n")
db <- build_custom_db(
    db_dir = file.path(opt$output_dir, "lincs_db"),
    overwrite = FALSE
)

# ============================================================================
# Load DE result files
# ============================================================================
cat("\nLoading DE results...\n")

de_files <- list.files(opt$de_results_dir, pattern = "\\.tsv$",
                       full.names = TRUE, recursive = TRUE)
cat(sprintf("  Found %d DE result files\n", length(de_files)))

if (length(de_files) == 0) {
    stop("No DE result files found in: ", opt$de_results_dir)
}

# ============================================================================
# Process each contrast
# ============================================================================
all_rankings <- list()

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

        # Get gene names
        gene_col <- colnames(de_df)[1]
        de_df$gene <- de_df[[gene_col]]

        # -----------------------------------------------------------------
        # 1. Extract query signature (top N up + top N down)
        # -----------------------------------------------------------------
        de_sig <- de_df %>%
            filter(!is.na(padj), !is.na(log2FoldChange))

        up_genes <- de_sig %>%
            filter(log2FoldChange > 0) %>%
            arrange(padj, desc(log2FoldChange)) %>%
            slice_head(n = opt$n_up) %>%
            pull(gene)

        down_genes <- de_sig %>%
            filter(log2FoldChange < 0) %>%
            arrange(padj, log2FoldChange) %>%
            slice_head(n = opt$n_down) %>%
            pull(gene)

        cat(sprintf("  Query signature: %d up, %d down genes\n",
                    length(up_genes), length(down_genes)))

        if (length(up_genes) < 10 || length(down_genes) < 10) {
            cat(sprintf("  WARNING: Too few query genes for %s, skipping\n",
                        contrast_name))
            next
        }

        # -----------------------------------------------------------------
        # 2. Run CMAP-style connectivity search
        # -----------------------------------------------------------------
        cat("  Running CMAP connectivity search against LINCS L1000...\n")

        # Create query signature
        qsig <- qSig(
            query = list(upset = up_genes, downset = down_genes),
            gess_method = "LINCS",
            refdb = db
        )

        # Run the GESS (Gene Expression Signature Search)
        gess_result <- gess_lincs(
            qSig = qsig,
            sortby = "NCS",
            tau = TRUE,
            workers = opt$threads
        )

        # Extract results
        drug_ranking <- result(gess_result)

        if (!is.null(drug_ranking) && nrow(drug_ranking) > 0) {
            cat(sprintf("  Found %d perturbagen results\n", nrow(drug_ranking)))

            # Sort by connectivity score (most negative = best reversal)
            drug_ranking <- drug_ranking %>%
                arrange(NCS)

            # Save full ranking
            write.table(drug_ranking,
                        file.path(ranking_dir,
                                  paste0(contrast_name, "_drug_ranking.tsv")),
                        sep = "\t", quote = FALSE, row.names = FALSE)

            # Save top candidates
            top_candidates <- head(drug_ranking, 100)
            write.table(top_candidates,
                        file.path(ranking_dir,
                                  paste0(contrast_name, "_top100_candidates.tsv")),
                        sep = "\t", quote = FALSE, row.names = FALSE)

            # Store for combined output
            top_candidates$contrast <- contrast_name
            all_rankings[[contrast_name]] <- top_candidates

            # -----------------------------------------------------------------
            # 3. Drug Set Enrichment Analysis (DSEA) by mechanism of action
            # -----------------------------------------------------------------
            cat("  Running Drug Set Enrichment Analysis (DSEA)...\n")

            tryCatch({
                dsea_result <- dsea_hyperG(
                    drugs = drug_ranking,
                    type = "MOA",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2
                )

                dsea_df <- result(dsea_result)

                if (!is.null(dsea_df) && nrow(dsea_df) > 0) {
                    write.table(dsea_df,
                                file.path(dsea_dir,
                                          paste0(contrast_name, "_dsea_moa.tsv")),
                                sep = "\t", quote = FALSE, row.names = FALSE)
                    cat(sprintf("  DSEA: %d significant MOA categories\n",
                                nrow(dsea_df)))

                    # DSEA dot plot
                    tryCatch({
                        plot_data <- head(dsea_df, 20)
                        if ("Description" %in% colnames(plot_data)) {
                            p_dsea <- ggplot(plot_data,
                                             aes(x = -log10(pvalue),
                                                 y = reorder(Description, -log10(pvalue)),
                                                 size = Count,
                                                 color = qvalue)) +
                                geom_point() +
                                scale_color_gradient(low = "red", high = "blue") +
                                theme_minimal() +
                                theme(axis.text.y = element_text(size = 8)) +
                                labs(
                                    title = paste("DSEA: MOA Enrichment -", contrast_name),
                                    x = "-log10(p-value)",
                                    y = "Mechanism of Action",
                                    size = "Gene Count",
                                    color = "q-value"
                                )
                            ggsave(file.path(plot_dir,
                                             paste0(contrast_name, "_dsea_dotplot.pdf")),
                                   p_dsea, width = 10, height = 8)
                        }
                    }, error = function(e) {
                        cat(sprintf("    WARNING: DSEA plot failed: %s\n",
                                    conditionMessage(e)))
                    })
                } else {
                    cat("  DSEA: No significant MOA enrichments found\n")
                }
            }, error = function(e) {
                cat(sprintf("  WARNING: DSEA failed: %s\n", conditionMessage(e)))
            })

            # -----------------------------------------------------------------
            # 4. Visualization: drug connectivity barplot
            # -----------------------------------------------------------------
            cat("  Generating drug ranking visualization...\n")

            tryCatch({
                # Determine perturbagen name column
                drug_name_col <- intersect(
                    c("pert", "pert_iname", "drug_name", "name"),
                    colnames(top_candidates)
                )
                if (length(drug_name_col) == 0) {
                    drug_name_col <- colnames(top_candidates)[1]
                } else {
                    drug_name_col <- drug_name_col[1]
                }

                plot_data <- head(top_candidates, 30)

                p_rank <- ggplot(plot_data,
                                 aes(x = reorder(!!sym(drug_name_col), -NCS),
                                     y = NCS)) +
                    geom_bar(stat = "identity",
                             fill = ifelse(plot_data$NCS < 0, "steelblue", "coral")) +
                    coord_flip() +
                    theme_minimal() +
                    theme(axis.text.y = element_text(size = 8)) +
                    labs(
                        title = paste("Top Drug Candidates -", contrast_name),
                        subtitle = "LINCS L1000 Normalized Connectivity Score (NCS)",
                        x = "Drug / Perturbagen",
                        y = "NCS (negative = reversal)"
                    ) +
                    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")

                ggsave(file.path(plot_dir,
                                 paste0(contrast_name, "_drug_ranking.pdf")),
                       p_rank, width = 10, height = 8)
            }, error = function(e) {
                cat(sprintf("    WARNING: Ranking plot failed: %s\n",
                            conditionMessage(e)))
            })

        } else {
            cat(sprintf("  WARNING: No results from GESS for %s\n", contrast_name))
        }

    }, error = function(e) {
        cat(sprintf("  ERROR processing %s: %s\n", contrast_name, conditionMessage(e)))
    })
}

# ============================================================================
# Combined summary across contrasts
# ============================================================================
cat("\n--- Combined Summary ---\n")

if (length(all_rankings) > 0) {
    combined <- bind_rows(all_rankings)
    write.table(combined,
                file.path(ranking_dir, "all_contrasts_top_candidates.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Combined results: %d drug candidates across %d contrasts\n",
                nrow(combined), length(all_rankings)))

    # Identify drugs appearing in multiple contrasts
    drug_name_col <- intersect(
        c("pert", "pert_iname", "drug_name", "name"),
        colnames(combined)
    )
    if (length(drug_name_col) > 0) {
        recurrent <- combined %>%
            count(!!sym(drug_name_col[1]), sort = TRUE) %>%
            filter(n > 1)

        if (nrow(recurrent) > 0) {
            write.table(recurrent,
                        file.path(ranking_dir, "recurrent_drug_candidates.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            cat(sprintf("  %d drugs appear in multiple contrasts\n",
                        nrow(recurrent)))
        }
    }
} else {
    cat("  WARNING: No drug rankings generated from any contrast\n")
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== signatureSearch Analysis Complete ===\n")
cat(sprintf("Output directory: %s\n", opt$output_dir))
cat("Key outputs:\n")
cat(sprintf("  - Drug rankings: %s\n", ranking_dir))
cat(sprintf("  - DSEA MOA enrichment: %s\n", dsea_dir))
cat(sprintf("  - Plots: %s\n", plot_dir))
