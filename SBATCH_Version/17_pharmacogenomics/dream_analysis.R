#!/usr/bin/env Rscript
# ============================================================================
# 17 - DREAM Druggability Analysis
# ============================================================================
# Implements network topology-based druggability scoring inspired by the DREAM
# (Drug Repurposing and Evaluation through Augmented Modules) approach:
#   1. Identify disease genes from DE results
#   2. Build co-expression network from normalized counts (Pearson correlation)
#   3. Compute network topology metrics (degree, betweenness, closeness)
#   4. Score gene druggability based on topological importance
#   5. Identify drug combination candidates from Louvain community modules
#
# Required packages: igraph, WGCNA, ggplot2, tidyverse, optparse, ggraph
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(igraph)
    library(tidyverse)
    library(ggplot2)
    library(ggraph)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
    make_option("--de_results_dir", type = "character",
                help = "Directory containing DE results (TSV files with log2FC and padj)"),
    make_option("--normalized_counts", type = "character",
                help = "Path to normalized count matrix (genes x samples TSV)"),
    make_option("--output_dir", type = "character",
                help = "Output directory"),
    make_option("--cor_threshold", type = "double", default = 0.7,
                help = "Correlation threshold for co-expression edges [default: %default]"),
    make_option("--padj_threshold", type = "double", default = 0.05,
                help = "Adjusted p-value threshold for DE genes [default: %default]"),
    make_option("--lfc_threshold", type = "double", default = 0.585,
                help = "Log2 fold-change threshold for DE genes [default: %default]"),
    make_option("--threads", type = "integer", default = 8,
                help = "Number of threads [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$de_results_dir) || is.null(opt$normalized_counts) ||
    is.null(opt$output_dir)) {
    stop("Required: --de_results_dir, --normalized_counts, --output_dir")
}

cat("=== DREAM Druggability Analysis ===\n")
cat("DE results dir:", opt$de_results_dir, "\n")
cat("Normalized counts:", opt$normalized_counts, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Correlation threshold:", opt$cor_threshold, "\n")
cat("padj threshold:", opt$padj_threshold, "\n")
cat("LFC threshold:", opt$lfc_threshold, "\n")
cat("Threads:", opt$threads, "\n\n")

# ============================================================================
# Create output directories
# ============================================================================
drug_score_dir <- file.path(opt$output_dir, "druggability_scores")
combo_dir <- file.path(opt$output_dir, "combination_candidates")
network_dir <- file.path(opt$output_dir, "networks")
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(drug_score_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Load normalized counts
# ============================================================================
cat("Loading normalized counts...\n")

counts <- read.table(opt$normalized_counts, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE)
cat(sprintf("  Count matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))

# ============================================================================
# Load DE results
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

        # Get gene name column
        gene_col <- colnames(de_df)[1]
        de_df$gene <- de_df[[gene_col]]

        # -----------------------------------------------------------------
        # 1. Identify disease genes
        # -----------------------------------------------------------------
        disease_genes <- de_df %>%
            filter(!is.na(padj), padj < opt$padj_threshold,
                   abs(log2FoldChange) > opt$lfc_threshold) %>%
            pull(gene)

        cat(sprintf("  Disease genes: %d\n", length(disease_genes)))

        if (length(disease_genes) < 20) {
            cat(sprintf("  WARNING: Too few disease genes for %s, skipping\n",
                        contrast_name))
            next
        }

        # Limit to top 500 disease genes for computational feasibility
        if (length(disease_genes) > 500) {
            top_de <- de_df %>%
                filter(gene %in% disease_genes) %>%
                arrange(padj) %>%
                slice_head(n = 500)
            disease_genes <- top_de$gene
            cat(sprintf("  Limited to top 500 disease genes\n"))
        }

        # Filter counts to disease genes present in expression data
        disease_genes <- disease_genes[disease_genes %in% rownames(counts)]
        cat(sprintf("  Disease genes in expression data: %d\n", length(disease_genes)))

        if (length(disease_genes) < 20) {
            cat(sprintf("  WARNING: Too few mapped genes, skipping\n"))
            next
        }

        # -----------------------------------------------------------------
        # 2. Build co-expression network
        # -----------------------------------------------------------------
        cat("  Building co-expression network...\n")

        # Subset expression matrix to disease genes
        expr_sub <- t(counts[disease_genes, , drop = FALSE])

        # Compute pairwise Pearson correlations
        cor_mat <- cor(expr_sub, method = "pearson", use = "pairwise.complete.obs")

        # Threshold to create adjacency
        adj_mat <- abs(cor_mat)
        adj_mat[adj_mat < opt$cor_threshold] <- 0
        diag(adj_mat) <- 0

        # Create igraph object
        g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected",
                                          weighted = TRUE)

        # Remove isolated vertices
        isolated <- which(degree(g) == 0)
        if (length(isolated) > 0) {
            g <- delete_vertices(g, isolated)
        }

        n_nodes <- vcount(g)
        n_edges <- ecount(g)
        cat(sprintf("  Network: %d nodes, %d edges\n", n_nodes, n_edges))

        if (n_nodes < 10) {
            cat(sprintf("  WARNING: Network too sparse for %s, skipping\n",
                        contrast_name))
            next
        }

        # -----------------------------------------------------------------
        # 3. Compute topology metrics
        # -----------------------------------------------------------------
        cat("  Computing network topology metrics...\n")

        node_degree <- degree(g)
        node_betweenness <- betweenness(g, directed = FALSE, normalized = TRUE)
        node_closeness <- closeness(g, mode = "all", normalized = TRUE)
        node_eigenvector <- tryCatch(
            eigen_centrality(g)$vector,
            error = function(e) rep(NA, n_nodes)
        )

        topology_df <- data.frame(
            gene = V(g)$name,
            degree = node_degree,
            betweenness = node_betweenness,
            closeness = node_closeness,
            eigenvector_centrality = node_eigenvector,
            stringsAsFactors = FALSE
        )

        # Add DE information
        topology_df <- topology_df %>%
            left_join(
                de_df %>% select(gene, log2FoldChange, padj),
                by = "gene"
            )

        # -----------------------------------------------------------------
        # 4. Compute druggability scores
        # -----------------------------------------------------------------
        cat("  Computing druggability scores...\n")

        # Normalize each metric to [0, 1]
        normalize_01 <- function(x) {
            x_clean <- x[!is.na(x)]
            if (length(x_clean) == 0 || max(x_clean) == min(x_clean)) return(rep(0, length(x)))
            (x - min(x_clean, na.rm = TRUE)) / (max(x_clean, na.rm = TRUE) - min(x_clean, na.rm = TRUE))
        }

        topology_df <- topology_df %>%
            mutate(
                degree_norm = normalize_01(degree),
                betweenness_norm = normalize_01(betweenness),
                closeness_norm = normalize_01(closeness),
                eigen_norm = normalize_01(eigenvector_centrality),
                lfc_norm = normalize_01(abs(log2FoldChange)),
                sig_norm = normalize_01(-log10(pmax(padj, 1e-300)))
            ) %>%
            mutate(
                # Composite druggability score: weighted combination of topology
                # and differential expression
                druggability_score = 0.20 * degree_norm +
                                     0.25 * betweenness_norm +
                                     0.15 * closeness_norm +
                                     0.15 * eigen_norm +
                                     0.15 * lfc_norm +
                                     0.10 * sig_norm
            ) %>%
            arrange(desc(druggability_score))

        # Save druggability scores
        score_file <- file.path(drug_score_dir,
                                paste0(contrast_name, "_druggability_scores.tsv"))
        write.table(topology_df, score_file,
                    sep = "\t", quote = FALSE, row.names = FALSE)
        cat(sprintf("  Saved druggability scores for %d genes\n", nrow(topology_df)))

        # -----------------------------------------------------------------
        # 5. Louvain community detection for combination candidates
        # -----------------------------------------------------------------
        cat("  Running Louvain community detection...\n")

        communities <- cluster_louvain(g)
        membership <- membership(communities)
        n_modules <- max(membership)
        modularity_score <- modularity(communities)

        cat(sprintf("  Detected %d modules (modularity: %.3f)\n",
                    n_modules, modularity_score))

        # Annotate genes with module membership
        topology_df$module <- membership[match(topology_df$gene, V(g)$name)]

        # Identify combination candidates: top gene per module
        combination_candidates <- topology_df %>%
            filter(!is.na(module)) %>%
            group_by(module) %>%
            arrange(desc(druggability_score)) %>%
            mutate(
                module_size = n(),
                rank_in_module = row_number()
            ) %>%
            ungroup()

        # Module summary: best target per module
        module_summary <- combination_candidates %>%
            group_by(module) %>%
            summarise(
                module_size = first(module_size),
                top_gene = gene[1],
                top_druggability = druggability_score[1],
                top_genes_3 = paste(head(gene, 3), collapse = "; "),
                mean_druggability = mean(druggability_score, na.rm = TRUE),
                mean_betweenness = mean(betweenness, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            arrange(desc(top_druggability))

        combo_file <- file.path(combo_dir,
                                paste0(contrast_name, "_combination_candidates.tsv"))
        write.table(combination_candidates, combo_file,
                    sep = "\t", quote = FALSE, row.names = FALSE)

        module_file <- file.path(combo_dir,
                                 paste0(contrast_name, "_module_summary.tsv"))
        write.table(module_summary, module_file,
                    sep = "\t", quote = FALSE, row.names = FALSE)

        cat(sprintf("  Saved combination candidates and %d module summaries\n",
                    nrow(module_summary)))

        # Save network as edge list
        edge_list <- as_data_frame(g, what = "edges")
        write.table(edge_list,
                    file.path(network_dir, paste0(contrast_name, "_edge_list.tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # -----------------------------------------------------------------
        # 6. Visualization
        # -----------------------------------------------------------------
        cat("  Generating visualizations...\n")

        # Druggability score barplot (top 30 genes)
        tryCatch({
            plot_data <- head(topology_df, 30)

            p_bar <- ggplot(plot_data,
                            aes(x = reorder(gene, druggability_score),
                                y = druggability_score,
                                fill = factor(module))) +
                geom_bar(stat = "identity") +
                coord_flip() +
                theme_minimal() +
                theme(axis.text.y = element_text(size = 8)) +
                labs(
                    title = paste("Top Druggable Genes -", contrast_name),
                    subtitle = "Composite score: topology + differential expression",
                    x = "Gene",
                    y = "Druggability Score",
                    fill = "Module"
                )
            ggsave(file.path(plot_dir,
                             paste0(contrast_name, "_druggability_barplot.pdf")),
                   p_bar, width = 10, height = 8)
        }, error = function(e) {
            cat(sprintf("    WARNING: Barplot failed: %s\n", conditionMessage(e)))
        })

        # Network visualization (top 100 nodes by druggability)
        tryCatch({
            top_nodes <- head(topology_df$gene, min(100, nrow(topology_df)))
            g_sub <- induced_subgraph(g, vids = which(V(g)$name %in% top_nodes))

            if (vcount(g_sub) > 5) {
                V(g_sub)$module <- as.character(membership[V(g_sub)$name])
                V(g_sub)$druggability <- topology_df$druggability_score[
                    match(V(g_sub)$name, topology_df$gene)
                ]

                p_net <- ggraph(g_sub, layout = "fr") +
                    geom_edge_link(alpha = 0.2, colour = "grey60") +
                    geom_node_point(aes(size = druggability, color = module),
                                   alpha = 0.8) +
                    geom_node_text(aes(label = name), size = 2, repel = TRUE) +
                    theme_void() +
                    labs(
                        title = paste("Co-Expression Network -", contrast_name),
                        subtitle = "Node size = druggability score, color = module",
                        size = "Druggability",
                        color = "Module"
                    )
                ggsave(file.path(plot_dir,
                                 paste0(contrast_name, "_network.pdf")),
                       p_net, width = 12, height = 10)
            }
        }, error = function(e) {
            cat(sprintf("    WARNING: Network plot failed: %s\n", conditionMessage(e)))
        })

    }, error = function(e) {
        cat(sprintf("  ERROR processing %s: %s\n", contrast_name, conditionMessage(e)))
    })
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== DREAM Druggability Analysis Complete ===\n")
cat(sprintf("Output directory: %s\n", opt$output_dir))
cat("Key outputs:\n")
cat(sprintf("  - Druggability scores: %s\n", drug_score_dir))
cat(sprintf("  - Combination candidates: %s\n", combo_dir))
cat(sprintf("  - Network edge lists: %s\n", network_dir))
cat(sprintf("  - Plots: %s\n", plot_dir))
