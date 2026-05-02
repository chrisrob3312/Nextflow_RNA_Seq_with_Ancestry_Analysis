#!/usr/bin/env Rscript
# ============================================================================
# Pathway Enrichment Analysis: ORA, GSEA, and GSVA
# ============================================================================
# For each DE contrast: runs clusterProfiler ORA on significant genes,
# fgsea GSEA on ranked gene lists, and GSVA for per-sample pathway scores.
# Queries MSigDB collections via msigdbr.
# ============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(fgsea)
  library(GSVA)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(optparse)
  library(BiocParallel)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
  make_option("--de_results_dir", type = "character",
              help = "Directory containing DE results (with deseq2/ and limma_voom/ subdirs)"),
  make_option("--norm_counts", type = "character",
              help = "Path to normalized count matrix for GSVA"),
  make_option("--metadata", type = "character",
              help = "Path to sample metadata"),
  make_option("--output_dir", type = "character",
              help = "Output directory"),
  make_option("--pathway_dbs", type = "character",
              default = "GO_BP,GO_MF,KEGG,REACTOME,HALLMARK,IMMUNESIGDB",
              help = "Comma-separated MSigDB collections to query"),
  make_option("--padj_threshold", type = "double", default = 0.05,
              help = "Adjusted p-value threshold for significant genes (ORA input)"),
  make_option("--lfc_threshold", type = "double", default = 0.585,
              help = "Log2FC threshold for significant genes"),
  make_option("--species", type = "character", default = "Homo sapiens",
              help = "Species for msigdbr"),
  make_option("--threads", type = "integer", default = 8,
              help = "Number of threads")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("de_results_dir", "norm_counts", "output_dir")) {
  if (is.null(opt[[arg]])) stop(sprintf("Required argument --%s is missing", arg))
}

register(MulticoreParam(workers = opt$threads))

cat("============================================================\n")
cat("Pathway Enrichment Analysis\n")
cat("============================================================\n")
cat(sprintf("DE results dir:  %s\n", opt$de_results_dir))
cat(sprintf("Norm counts:     %s\n", opt$norm_counts))
cat(sprintf("Metadata:        %s\n", opt$metadata))
cat(sprintf("Output dir:      %s\n", opt$output_dir))
cat(sprintf("Pathway DBs:     %s\n", opt$pathway_dbs))
cat(sprintf("padj threshold:  %s\n", opt$padj_threshold))
cat(sprintf("LFC threshold:   %s\n", opt$lfc_threshold))
cat(sprintf("Species:         %s\n", opt$species))
cat(sprintf("Threads:         %d\n", opt$threads))
cat("============================================================\n\n")

# ============================================================================
# [1] Load MSigDB gene sets
# ============================================================================
cat("[1] Loading MSigDB gene sets...\n")

pathway_dbs <- trimws(unlist(strsplit(opt$pathway_dbs, ",")))
cat(sprintf("    Requested collections: %s\n", paste(pathway_dbs, collapse = ", ")))

# Map user-friendly names to msigdbr categories/subcategories
db_mapping <- list(
  GO_BP = list(category = "C5", subcategory = "GO:BP"),
  GO_MF = list(category = "C5", subcategory = "GO:MF"),
  GO_CC = list(category = "C5", subcategory = "GO:CC"),
  KEGG = list(category = "C2", subcategory = "CP:KEGG"),
  REACTOME = list(category = "C2", subcategory = "CP:REACTOME"),
  HALLMARK = list(category = "H", subcategory = NA),
  IMMUNESIGDB = list(category = "C7", subcategory = "IMMUNESIGDB")
)

# Fetch gene sets
all_gene_sets <- list()
gene_set_lists <- list()  # For fgsea (named list of gene vectors)

for (db_name in pathway_dbs) {
  if (!db_name %in% names(db_mapping)) {
    cat(sprintf("    WARNING: Unknown database '%s', skipping.\n", db_name))
    next
  }

  db_info <- db_mapping[[db_name]]
  cat(sprintf("    Fetching %s...\n", db_name))

  if (is.na(db_info$subcategory)) {
    msig_df <- msigdbr(species = opt$species, category = db_info$category)
  } else {
    msig_df <- msigdbr(species = opt$species, category = db_info$category,
                       subcategory = db_info$subcategory)
  }

  if (nrow(msig_df) == 0) {
    cat(sprintf("    WARNING: No gene sets found for %s\n", db_name))
    next
  }

  # Store as data frame for clusterProfiler
  all_gene_sets[[db_name]] <- msig_df[, c("gs_name", "entrez_gene", "gene_symbol")]

  # Create named list for fgsea (gene symbols)
  gs_list <- split(msig_df$gene_symbol, msig_df$gs_name)
  gene_set_lists[[db_name]] <- gs_list

  cat(sprintf("      %d gene sets, %d unique genes\n",
              length(unique(msig_df$gs_name)), length(unique(msig_df$gene_symbol))))
}

# ============================================================================
# [2] Load normalized counts for GSVA
# ============================================================================
cat("\n[2] Loading normalized counts for GSVA...\n")

norm_counts <- read.delim(opt$norm_counts, row.names = 1, check.names = FALSE)
cat(sprintf("    Matrix: %d genes x %d samples\n", nrow(norm_counts), ncol(norm_counts)))

# Load metadata if provided
if (!is.null(opt$metadata) && file.exists(opt$metadata)) {
  metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE)
  if (!"sample_id" %in% colnames(metadata)) {
    colnames(metadata)[1] <- "sample_id"
  }
  rownames(metadata) <- metadata$sample_id
} else {
  metadata <- data.frame(sample_id = colnames(norm_counts), row.names = colnames(norm_counts))
}

# ============================================================================
# [3] Find DE result files
# ============================================================================
cat("\n[3] Finding DE result files...\n")

# Look for DESeq2 results primarily
de_files <- list.files(file.path(opt$de_results_dir, "deseq2"),
                       pattern = "\\.tsv$", full.names = TRUE)

if (length(de_files) == 0) {
  # Fallback to limma_voom
  de_files <- list.files(file.path(opt$de_results_dir, "limma_voom"),
                         pattern = "\\.tsv$", full.names = TRUE)
}

if (length(de_files) == 0) {
  stop("ERROR: No DE result files found in ", opt$de_results_dir)
}

cat(sprintf("    Found %d DE result files\n", length(de_files)))

# ============================================================================
# Helper Functions
# ============================================================================

run_ora <- function(sig_genes, universe_genes, gene_set_df, db_name, contrast_name,
                    direction, output_dir) {
  # Over-Representation Analysis using clusterProfiler enricher()
  if (length(sig_genes) < 5) {
    cat(sprintf("        Skipping ORA (%s, %s): <5 significant genes\n", db_name, direction))
    return(NULL)
  }

  # Use gene symbols for enrichment
  term2gene <- gene_set_df[, c("gs_name", "gene_symbol")]
  colnames(term2gene) <- c("term", "gene")

  result <- tryCatch({
    enricher(
      gene = sig_genes,
      universe = universe_genes,
      TERM2GENE = term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1,
      minGSSize = 10,
      maxGSSize = 500
    )
  }, error = function(e) {
    cat(sprintf("        ORA error (%s): %s\n", db_name, conditionMessage(e)))
    return(NULL)
  })

  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    cat(sprintf("        No significant ORA results for %s (%s)\n", db_name, direction))
    return(NULL)
  }

  # Save results
  res_df <- as.data.frame(result)
  out_file <- file.path(output_dir, "ora",
                        sprintf("ora_%s_%s_%s.tsv",
                                gsub("[^A-Za-z0-9_]", "_", contrast_name),
                                db_name, direction))
  write.table(res_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

  return(result)
}

run_gsea_fgsea <- function(ranked_genes, gene_set_list, db_name, contrast_name, output_dir) {
  # GSEA using fgsea
  if (length(ranked_genes) < 100) {
    cat(sprintf("        Skipping GSEA (%s): <100 ranked genes\n", db_name))
    return(NULL)
  }

  result <- tryCatch({
    fgsea(
      pathways = gene_set_list,
      stats = ranked_genes,
      minSize = 15,
      maxSize = 500,
      nPermSimple = 10000
    )
  }, error = function(e) {
    cat(sprintf("        GSEA error (%s): %s\n", db_name, conditionMessage(e)))
    return(NULL)
  })

  if (is.null(result) || nrow(result) == 0) {
    cat(sprintf("        No GSEA results for %s\n", db_name))
    return(NULL)
  }

  # Sort by NES

  result <- result[order(result$padj), ]

  # Save results (convert leadingEdge list to string for TSV)
  res_out <- result
  res_out$leadingEdge <- sapply(res_out$leadingEdge, paste, collapse = ";")
  out_file <- file.path(output_dir, "gsea",
                        sprintf("gsea_%s_%s.tsv",
                                gsub("[^A-Za-z0-9_]", "_", contrast_name), db_name))
  write.table(res_out, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

  return(result)
}

generate_ora_dotplot <- function(ora_result, db_name, contrast_name, direction, output_dir) {
  if (is.null(ora_result) || nrow(as.data.frame(ora_result)) == 0) return(invisible(NULL))

  res_df <- as.data.frame(ora_result)
  n_show <- min(20, nrow(res_df))
  res_df <- head(res_df[order(res_df$p.adjust), ], n_show)

  # Compute gene ratio as numeric
  res_df$GeneRatioNumeric <- sapply(res_df$GeneRatio, function(x) {
    parts <- as.numeric(unlist(strsplit(x, "/")))
    parts[1] / parts[2]
  })

  p <- ggplot(res_df, aes(x = GeneRatioNumeric, y = reorder(Description, GeneRatioNumeric),
                           size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = sprintf("ORA: %s (%s, %s)", contrast_name, db_name, direction),
         x = "Gene Ratio", y = "") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))

  filename <- file.path(output_dir, "plots",
                        sprintf("dotplot_ora_%s_%s_%s.pdf",
                                gsub("[^A-Za-z0-9_]", "_", contrast_name), db_name, direction))
  ggsave(filename, p, width = 10, height = max(4, n_show * 0.35))
}

generate_gsea_plot <- function(gsea_result, db_name, contrast_name, output_dir) {
  if (is.null(gsea_result) || nrow(gsea_result) == 0) return(invisible(NULL))

  # Barplot of top enriched/depleted pathways by NES
  sig_results <- gsea_result[gsea_result$padj < 0.05, ]
  if (nrow(sig_results) == 0) {
    sig_results <- head(gsea_result[order(gsea_result$pval), ], 20)
  }

  # Top positive and negative NES
  top_pos <- head(sig_results[order(-sig_results$NES), ], 10)
  top_neg <- head(sig_results[order(sig_results$NES), ], 10)
  plot_df <- rbind(top_pos, top_neg)
  plot_df <- plot_df[!duplicated(plot_df$pathway), ]

  if (nrow(plot_df) < 2) return(invisible(NULL))

  plot_df$direction <- ifelse(plot_df$NES > 0, "Activated", "Suppressed")

  p <- ggplot(plot_df, aes(x = NES, y = reorder(pathway, NES), fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Activated" = "firebrick", "Suppressed" = "steelblue")) +
    labs(title = sprintf("GSEA: %s (%s)", contrast_name, db_name),
         x = "Normalized Enrichment Score", y = "") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 7))

  filename <- file.path(output_dir, "plots",
                        sprintf("barplot_gsea_%s_%s.pdf",
                                gsub("[^A-Za-z0-9_]", "_", contrast_name), db_name))
  ggsave(filename, p, width = 10, height = max(4, nrow(plot_df) * 0.35))
}

# ============================================================================
# [4] Process each contrast
# ============================================================================
cat("\n[4] Running pathway enrichment for each contrast...\n")
cat("============================================================\n")

all_ora_results <- list()
all_gsea_results <- list()

for (de_file in de_files) {
  contrast_name <- gsub("\\.tsv$", "", basename(de_file))
  cat(sprintf("\n--- Processing: %s ---\n", contrast_name))

  # Load DE results
  de_res <- tryCatch({
    read.delim(de_file, row.names = 1, check.names = FALSE)
  }, error = function(e) {
    cat(sprintf("    ERROR reading file: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(de_res)) next

  # Standardize column names
  if ("log2FoldChange" %in% colnames(de_res)) {
    lfc_col <- "log2FoldChange"
  } else if ("logFC" %in% colnames(de_res)) {
    lfc_col <- "logFC"
  } else {
    cat("    WARNING: Cannot find LFC column, skipping.\n")
    next
  }

  if ("padj" %in% colnames(de_res)) {
    padj_col <- "padj"
  } else if ("adj.P.Val" %in% colnames(de_res)) {
    padj_col <- "adj.P.Val"
  } else {
    cat("    WARNING: Cannot find padj column, skipping.\n")
    next
  }

  # Define gene sets
  all_genes <- rownames(de_res)
  sig_up <- rownames(de_res)[!is.na(de_res[[padj_col]]) &
                               de_res[[padj_col]] < opt$padj_threshold &
                               de_res[[lfc_col]] > opt$lfc_threshold]
  sig_down <- rownames(de_res)[!is.na(de_res[[padj_col]]) &
                                 de_res[[padj_col]] < opt$padj_threshold &
                                 de_res[[lfc_col]] < -opt$lfc_threshold]
  sig_all <- c(sig_up, sig_down)

  cat(sprintf("    Genes: %d total, %d up, %d down\n",
              length(all_genes), length(sig_up), length(sig_down)))

  # Create ranked gene list for GSEA (by -log10(pval) * sign(LFC))
  pval_col <- ifelse("pvalue" %in% colnames(de_res), "pvalue",
                     ifelse("P.Value" %in% colnames(de_res), "P.Value", padj_col))

  ranked_genes <- de_res[[lfc_col]]
  names(ranked_genes) <- rownames(de_res)
  # Use stat if available (more informative), otherwise use signed -log10(pval)
  if ("stat" %in% colnames(de_res)) {
    ranked_genes <- de_res$stat
    names(ranked_genes) <- rownames(de_res)
  } else if ("t_statistic" %in% colnames(de_res)) {
    ranked_genes <- de_res$t_statistic
    names(ranked_genes) <- rownames(de_res)
  }
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  # Run enrichment for each database

  for (db_name in names(all_gene_sets)) {
    cat(sprintf("    [%s]\n", db_name))
    gene_set_df <- all_gene_sets[[db_name]]
    gs_list <- gene_set_lists[[db_name]]

    # --- ORA: up-regulated genes ---
    ora_up <- run_ora(sig_up, all_genes, gene_set_df, db_name, contrast_name, "up", opt$output_dir)
    if (!is.null(ora_up)) {
      generate_ora_dotplot(ora_up, db_name, contrast_name, "up", opt$output_dir)
      all_ora_results[[paste(contrast_name, db_name, "up", sep = "_")]] <- ora_up
    }

    # --- ORA: down-regulated genes ---
    ora_down <- run_ora(sig_down, all_genes, gene_set_df, db_name, contrast_name, "down", opt$output_dir)
    if (!is.null(ora_down)) {
      generate_ora_dotplot(ora_down, db_name, contrast_name, "down", opt$output_dir)
      all_ora_results[[paste(contrast_name, db_name, "down", sep = "_")]] <- ora_down
    }

    # --- ORA: all significant genes ---
    ora_all <- run_ora(sig_all, all_genes, gene_set_df, db_name, contrast_name, "all", opt$output_dir)
    if (!is.null(ora_all)) {
      generate_ora_dotplot(ora_all, db_name, contrast_name, "all", opt$output_dir)
    }

    # --- GSEA (fgsea) ---
    gsea_res <- run_gsea_fgsea(ranked_genes, gs_list, db_name, contrast_name, opt$output_dir)
    if (!is.null(gsea_res)) {
      generate_gsea_plot(gsea_res, db_name, contrast_name, opt$output_dir)
      all_gsea_results[[paste(contrast_name, db_name, sep = "_")]] <- gsea_res
    }
  }
}

# ============================================================================
# [5] GSVA - Per-sample pathway scores
# ============================================================================
cat("\n\n[5] Running GSVA for per-sample pathway scores...\n")

# Use HALLMARK gene sets for GSVA (manageable size)
gsva_db <- "HALLMARK"
if (gsva_db %in% names(gene_set_lists)) {
  gs_for_gsva <- gene_set_lists[[gsva_db]]
  cat(sprintf("    Using %d %s gene sets for GSVA\n", length(gs_for_gsva), gsva_db))

  # Filter gene sets to genes present in the expression matrix
  gs_for_gsva <- lapply(gs_for_gsva, function(genes) {
    intersect(genes, rownames(norm_counts))
  })
  gs_for_gsva <- gs_for_gsva[sapply(gs_for_gsva, length) >= 10]
  cat(sprintf("    Gene sets with >= 10 genes in matrix: %d\n", length(gs_for_gsva)))

  if (length(gs_for_gsva) >= 2) {
    # Run GSVA
    expr_matrix <- as.matrix(norm_counts)

    gsva_scores <- tryCatch({
      gsva(
        expr_matrix,
        gs_for_gsva,
        method = "gsva",
        kcdf = "Gaussian",
        parallel.sz = opt$threads,
        verbose = FALSE
      )
    }, error = function(e) {
      cat(sprintf("    GSVA error: %s\n", conditionMessage(e)))
      # Try with gsvaParam for newer GSVA versions
      tryCatch({
        param <- gsvaParam(expr_matrix, gs_for_gsva, kcdf = "Gaussian")
        gsva(param, BPPARAM = MulticoreParam(workers = opt$threads))
      }, error = function(e2) {
        cat(sprintf("    GSVA retry error: %s\n", conditionMessage(e2)))
        return(NULL)
      })
    })

    if (!is.null(gsva_scores)) {
      cat(sprintf("    GSVA computed: %d pathways x %d samples\n",
                  nrow(gsva_scores), ncol(gsva_scores)))

      # Save GSVA scores
      gsva_df <- as.data.frame(gsva_scores)
      gsva_df$pathway <- rownames(gsva_df)
      gsva_df <- gsva_df[, c("pathway", setdiff(colnames(gsva_df), "pathway"))]
      write.table(gsva_df,
                  file = file.path(opt$output_dir, "gsva", "gsva_hallmark_scores.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      # Generate GSVA heatmap
      cat("    Generating GSVA heatmap...\n")

      # Prepare annotation
      anno_cols <- intersect(c("condition", "group", "predicted_ancestry", "subtype"),
                             colnames(metadata))
      common_gsva_samples <- intersect(colnames(gsva_scores), rownames(metadata))

      if (length(anno_cols) > 0 && length(common_gsva_samples) > 0) {
        annotation_col <- metadata[common_gsva_samples, anno_cols, drop = FALSE]
      } else {
        annotation_col <- NA
      }

      # Heatmap of GSVA scores
      pdf(file.path(opt$output_dir, "plots", "gsva_hallmark_heatmap.pdf"),
          width = max(10, ncol(gsva_scores) * 0.15),
          height = max(8, nrow(gsva_scores) * 0.25))
      pheatmap(
        gsva_scores[, common_gsva_samples, drop = FALSE],
        annotation_col = annotation_col,
        show_colnames = (ncol(gsva_scores) <= 50),
        clustering_method = "ward.D2",
        color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
        main = sprintf("GSVA Scores: %s Pathways", gsva_db),
        fontsize_row = 7
      )
      dev.off()
      cat("    Saved: gsva_hallmark_heatmap.pdf\n")
    }
  } else {
    cat("    Not enough gene sets for GSVA, skipping.\n")
  }
} else {
  cat("    HALLMARK gene sets not loaded, skipping GSVA.\n")
}

# Also run GSVA for other collections if they are small enough
for (db_name in setdiff(names(gene_set_lists), gsva_db)) {
  gs_list <- gene_set_lists[[db_name]]
  # Only run GSVA for smaller collections (< 300 gene sets)
  if (length(gs_list) > 300) {
    cat(sprintf("    Skipping GSVA for %s (%d gene sets, too large)\n", db_name, length(gs_list)))
    next
  }

  cat(sprintf("    Running GSVA for %s (%d gene sets)...\n", db_name, length(gs_list)))

  gs_filtered <- lapply(gs_list, function(genes) intersect(genes, rownames(norm_counts)))
  gs_filtered <- gs_filtered[sapply(gs_filtered, length) >= 10]

  if (length(gs_filtered) < 2) next

  gsva_res <- tryCatch({
    gsva(as.matrix(norm_counts), gs_filtered, method = "gsva",
         kcdf = "Gaussian", parallel.sz = opt$threads, verbose = FALSE)
  }, error = function(e) {
    tryCatch({
      param <- gsvaParam(as.matrix(norm_counts), gs_filtered, kcdf = "Gaussian")
      gsva(param, BPPARAM = MulticoreParam(workers = opt$threads))
    }, error = function(e2) NULL)
  })

  if (!is.null(gsva_res)) {
    gsva_out <- as.data.frame(gsva_res)
    gsva_out$pathway <- rownames(gsva_out)
    gsva_out <- gsva_out[, c("pathway", setdiff(colnames(gsva_out), "pathway"))]
    write.table(gsva_out,
                file = file.path(opt$output_dir, "gsva",
                                 sprintf("gsva_%s_scores.tsv", tolower(db_name))),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("    Saved GSVA scores for %s\n", db_name))
  }
}

# ============================================================================
# [6] Save R objects and summary
# ============================================================================
cat("\n[6] Saving results...\n")

save(all_ora_results, all_gsea_results,
     file = file.path(opt$output_dir, "rdata", "pathway_results.RData"))
cat("    Saved: pathway_results.RData\n")

# Generate summary
summary_rows <- list()
for (name in names(all_gsea_results)) {
  res <- all_gsea_results[[name]]
  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  n_pos <- sum(res$padj < 0.05 & res$NES > 0, na.rm = TRUE)
  n_neg <- sum(res$padj < 0.05 & res$NES < 0, na.rm = TRUE)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    analysis = name,
    method = "GSEA",
    n_tested = nrow(res),
    n_significant = n_sig,
    n_activated = n_pos,
    n_suppressed = n_neg,
    stringsAsFactors = FALSE
  )
}

for (name in names(all_ora_results)) {
  res_df <- as.data.frame(all_ora_results[[name]])
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    analysis = name,
    method = "ORA",
    n_tested = nrow(res_df),
    n_significant = sum(res_df$p.adjust < 0.05, na.rm = TRUE),
    n_activated = NA,
    n_suppressed = NA,
    stringsAsFactors = FALSE
  )
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  write.table(summary_df,
              file = file.path(opt$output_dir, "pathway_enrichment_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("    Saved: pathway_enrichment_summary.tsv\n")
  cat("\n    Summary:\n")
  print(summary_df)
}

cat("\n============================================================\n")
cat("Pathway enrichment analysis complete!\n")
cat("============================================================\n")

# Session info
writeLines(capture.output(sessionInfo()),
           file.path(opt$output_dir, "session_info.txt"))
