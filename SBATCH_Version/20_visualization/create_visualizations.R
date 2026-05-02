#!/usr/bin/env Rscript
# ============================================================================
# INTEGRATED PIPELINE VISUALIZATIONS
# ============================================================================
# Generates summary figures from all pipeline outputs and creates an HTML
# index with thumbnails for easy browsing.
#
# Usage:
#   Rscript create_visualizations.R \
#     --results_dir <path> \
#     --output_dir <path> \
#     --threads <int>
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(htmltools)
})

# Try to load optional packages
has_complex_heatmap <- suppressWarnings(requireNamespace("ComplexHeatmap", quietly = TRUE))
if (has_complex_heatmap) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
}

# ============================================================================
# PARSE ARGUMENTS
# ============================================================================

option_list <- list(
  make_option("--results_dir", type = "character", default = NULL,
              help = "Root results directory containing all pipeline outputs"),
  make_option("--output_dir", type = "character", default = "visualization_output",
              help = "Output directory for figures and HTML"),
  make_option("--threads", type = "integer", default = 4,
              help = "Number of threads [default: 4]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Setup directories
figures_dir <- file.path(opt$output_dir, "figures")
html_dir <- file.path(opt$output_dir, "html")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)

# Track generated figures for HTML index
figure_registry <- data.frame(
  filename = character(),
  title = character(),
  category = character(),
  description = character(),
  stringsAsFactors = FALSE
)

register_figure <- function(filename, title, category, description) {
  figure_registry <<- rbind(figure_registry, data.frame(
    filename = filename,
    title = title,
    category = category,
    description = description,
    stringsAsFactors = FALSE
  ))
}

# ============================================================================
# HELPER: SAFELY LOAD FILES
# ============================================================================

safe_fread <- function(path, ...) {
  if (file.exists(path)) {
    tryCatch(fread(path, ...), error = function(e) {
      cat(sprintf("  WARNING: Failed to read %s: %s\n", path, e$message))
      return(NULL)
    })
  } else {
    NULL
  }
}

# ============================================================================
# 1. DIFFERENTIAL EXPRESSION - COMBINED VOLCANO OVERVIEW
# ============================================================================

cat("=== Generating DE Figures ===\n")

create_de_volcano <- function() {
  de_dir <- file.path(opt$results_dir, "07_differential_expression", "deseq2")
  if (!dir.exists(de_dir)) {
    de_dir <- file.path(opt$results_dir, "differential_expression", "deseq2")
  }
  if (!dir.exists(de_dir)) {
    cat("  DE results not found, skipping volcano plot.\n")
    return(NULL)
  }

  result_files <- list.files(de_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(result_files) == 0) return(NULL)

  all_de <- list()
  for (f in result_files) {
    res <- safe_fread(f)
    if (is.null(res)) next

    contrast_name <- gsub("_deseq2\\.tsv$|\\.tsv$", "", basename(f))

    # Identify columns
    lfc_col <- grep("log2FoldChange|logFC", names(res), value = TRUE, ignore.case = TRUE)[1]
    padj_col <- grep("padj|adj.P.Val|FDR", names(res), value = TRUE, ignore.case = TRUE)[1]

    if (is.na(lfc_col) || is.na(padj_col)) next

    df <- data.frame(
      lfc = res[[lfc_col]],
      padj = res[[padj_col]],
      contrast = contrast_name,
      stringsAsFactors = FALSE
    )
    df <- df[complete.cases(df), ]
    all_de[[contrast_name]] <- df
  }

  if (length(all_de) == 0) return(NULL)

  combined <- do.call(rbind, all_de)
  combined$sig <- combined$padj < 0.05 & abs(combined$lfc) > 0.585
  combined$direction <- ifelse(combined$lfc > 0, "Up", "Down")
  combined$direction[!combined$sig] <- "NS"

  # Faceted volcano
  p <- ggplot(combined, aes(x = lfc, y = -log10(padj), color = direction)) +
    geom_point(size = 0.5, alpha = 0.4) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70"),
                      name = "Direction") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "grey40", linewidth = 0.3) +
    facet_wrap(~contrast, scales = "free_y") +
    labs(
      title = "Differential Expression Overview",
      subtitle = "Combined volcano plots across all contrasts",
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~adjusted~p)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      legend.position = "bottom"
    )

  n_facets <- length(unique(combined$contrast))
  fig_width <- min(16, max(8, 4 * ceiling(sqrt(n_facets))))
  fig_height <- min(12, max(6, 4 * ceiling(n_facets / ceiling(sqrt(n_facets)))))

  outfile <- file.path(figures_dir, "de_volcano_overview.png")
  ggsave(outfile, p, width = fig_width, height = fig_height, dpi = 150)
  register_figure("de_volcano_overview.png", "DE Volcano Overview",
                  "Differential Expression",
                  sprintf("Combined volcano plots for %d contrasts", n_facets))
  cat(sprintf("  Generated: de_volcano_overview.png (%d contrasts)\n", n_facets))
}

create_de_volcano()

# ============================================================================
# 2. ANCESTRY - PCA AND PROPORTION BARS
# ============================================================================

cat("=== Generating Ancestry Figures ===\n")

create_ancestry_figures <- function() {
  anc_dir <- file.path(opt$results_dir, "04_ancestry")
  if (!dir.exists(anc_dir)) {
    anc_dir <- file.path(opt$results_dir, "ancestry")
  }
  if (!dir.exists(anc_dir)) {
    cat("  Ancestry results not found, skipping.\n")
    return(NULL)
  }

  # Look for ancestry results files
  anc_files <- list.files(anc_dir, pattern = "ancestry|graf|somalier", full.names = TRUE, recursive = TRUE)
  prop_file <- list.files(anc_dir, pattern = "proportions|ancestry_results", full.names = TRUE)[1]

  if (is.na(prop_file) || !file.exists(prop_file)) {
    cat("  No ancestry proportions file found, skipping.\n")
    return(NULL)
  }

  anc_data <- safe_fread(prop_file)
  if (is.null(anc_data)) return(NULL)

  # Try to identify ancestry/population column
  pop_col <- grep("population|ancestry|category|predicted|graf", names(anc_data),
                  value = TRUE, ignore.case = TRUE)[1]
  sample_col <- names(anc_data)[1]

  if (is.na(pop_col)) {
    cat("  Cannot identify ancestry column, skipping.\n")
    return(NULL)
  }

  # Look for PCA coordinates
  pc1_col <- grep("^PC1$|^pc1$|^pca1$", names(anc_data), value = TRUE, ignore.case = TRUE)[1]
  pc2_col <- grep("^PC2$|^pc2$|^pca2$", names(anc_data), value = TRUE, ignore.case = TRUE)[1]

  # PCA plot colored by ancestry

  if (!is.na(pc1_col) && !is.na(pc2_col)) {
    p_pca <- ggplot(anc_data, aes_string(x = pc1_col, y = pc2_col, color = pop_col)) +
      geom_point(size = 2.5, alpha = 0.7) +
      labs(
        title = "Ancestry PCA",
        subtitle = "Samples colored by GRAF-predicted ancestry",
        x = "PC1",
        y = "PC2",
        color = "Ancestry"
      ) +
      theme_minimal() +
      theme(legend.position = "right")

    outfile <- file.path(figures_dir, "ancestry_pca.png")
    ggsave(outfile, p_pca, width = 9, height = 7, dpi = 150)
    register_figure("ancestry_pca.png", "Ancestry PCA",
                    "Ancestry", "PCA colored by GRAF-predicted ancestry category")
    cat("  Generated: ancestry_pca.png\n")
  }

  # Stacked proportion bar chart
  pop_counts <- as.data.frame(table(anc_data[[pop_col]]))
  names(pop_counts) <- c("Ancestry", "Count")
  pop_counts$Proportion <- pop_counts$Count / sum(pop_counts$Count)

  p_bar <- ggplot(pop_counts, aes(x = "Cohort", y = Proportion, fill = Ancestry)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", Count, Proportion * 100)),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Ancestry Composition",
      subtitle = sprintf("N = %d samples", sum(pop_counts$Count)),
      y = "Proportion",
      x = ""
    ) +
    theme_minimal() +
    theme(axis.text.x = element_blank())

  outfile <- file.path(figures_dir, "ancestry_proportions.png")
  ggsave(outfile, p_bar, width = 6, height = 7, dpi = 150)
  register_figure("ancestry_proportions.png", "Ancestry Proportions",
                  "Ancestry", "Stacked bar chart of ancestry composition")
  cat("  Generated: ancestry_proportions.png\n")
}

create_ancestry_figures()

# ============================================================================
# 3. IMMUNE - CELL COMPOSITION SUMMARY
# ============================================================================

cat("=== Generating Immune Figures ===\n")

create_immune_figures <- function() {
  immune_dir <- file.path(opt$results_dir, "13_immune_analysis")
  if (!dir.exists(immune_dir)) {
    immune_dir <- file.path(opt$results_dir, "immune_analysis")
  }
  if (!dir.exists(immune_dir)) {
    cat("  Immune analysis results not found, skipping.\n")
    return(NULL)
  }

  # Look for immune deconvolution results
  immune_files <- list.files(immune_dir, pattern = "cibersort|xcell|epic|mcp|quantiseq|timer",
                            full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

  if (length(immune_files) == 0) {
    # Try generic names
    immune_files <- list.files(immune_dir, pattern = "\\.tsv$|\\.csv$", full.names = TRUE)
  }

  if (length(immune_files) == 0) {
    cat("  No immune deconvolution files found, skipping.\n")
    return(NULL)
  }

  all_immune <- list()
  for (f in immune_files) {
    method_name <- gsub("\\.tsv$|\\.csv$|_results$", "", basename(f))
    dat <- safe_fread(f)
    if (is.null(dat) || ncol(dat) < 2) next

    # Assume first column is sample, rest are cell types
    sample_col <- names(dat)[1]
    cell_cols <- names(dat)[-1]

    # Melt to long format
    dat_long <- melt(dat, id.vars = sample_col, variable.name = "cell_type",
                     value.name = "proportion")
    dat_long$method <- method_name
    names(dat_long)[1] <- "sample"
    all_immune[[method_name]] <- dat_long
  }

  if (length(all_immune) == 0) return(NULL)

  combined <- do.call(rbind, all_immune)

  # Summary: mean proportion by cell type and method
  summary_df <- combined[, .(mean_prop = mean(proportion, na.rm = TRUE),
                            sd_prop = sd(proportion, na.rm = TRUE)),
                        by = .(cell_type, method)]

  # Top cell types by mean proportion
  top_cells <- summary_df[, .(max_prop = max(mean_prop, na.rm = TRUE)), by = cell_type]
  top_cells <- top_cells[order(-max_prop)][1:min(15, nrow(top_cells))]

  summary_plot <- summary_df[cell_type %in% top_cells$cell_type]

  p <- ggplot(summary_plot, aes(x = reorder(cell_type, mean_prop), y = mean_prop, fill = method)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, mean_prop - sd_prop),
                      ymax = mean_prop + sd_prop),
                  position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.3) +
    coord_flip() +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Immune Cell Composition Summary",
      subtitle = sprintf("Top %d cell types across %d method(s)", nrow(top_cells), length(all_immune)),
      x = "Cell Type",
      y = "Mean Proportion",
      fill = "Method"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")

  outfile <- file.path(figures_dir, "immune_composition.png")
  ggsave(outfile, p, width = 10, height = 8, dpi = 150)
  register_figure("immune_composition.png", "Immune Cell Composition",
                  "Immune Analysis", "Summary of immune cell proportions across methods")
  cat("  Generated: immune_composition.png\n")
}

create_immune_figures()

# ============================================================================
# 4. WGCNA - MODULE-TRAIT CORRELATION SUMMARY
# ============================================================================

cat("=== Generating WGCNA Figures ===\n")

create_wgcna_figures <- function() {
  wgcna_dir <- file.path(opt$results_dir, "11_wgcna")
  if (!dir.exists(wgcna_dir)) {
    wgcna_dir <- file.path(opt$results_dir, "wgcna")
  }
  if (!dir.exists(wgcna_dir)) {
    cat("  WGCNA results not found, skipping.\n")
    return(NULL)
  }

  # Look for module-trait correlation file
  mt_file <- list.files(wgcna_dir, pattern = "module.trait|moduleTrait|module_trait",
                       full.names = TRUE, recursive = TRUE, ignore.case = TRUE)[1]

  if (is.na(mt_file) || !file.exists(mt_file)) {
    # Try RData
    rdata_file <- list.files(wgcna_dir, pattern = "\\.RData$|\\.rds$",
                            full.names = TRUE, recursive = TRUE)[1]
    if (!is.na(rdata_file)) {
      cat("  Found WGCNA RData but no module-trait TSV. Skipping heatmap.\n")
    } else {
      cat("  No WGCNA module-trait results found, skipping.\n")
    }
    return(NULL)
  }

  mt_data <- safe_fread(mt_file)
  if (is.null(mt_data)) return(NULL)

  # Expect: rows = modules, columns = traits, values = correlations
  # First column might be module names
  if (!is.numeric(mt_data[[1]])) {
    module_names <- mt_data[[1]]
    mt_matrix <- as.matrix(mt_data[, -1, with = FALSE])
    rownames(mt_matrix) <- module_names
  } else {
    mt_matrix <- as.matrix(mt_data)
  }

  # Use ComplexHeatmap if available, otherwise ggplot
  if (has_complex_heatmap) {
    col_fun <- colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C"))

    outfile <- file.path(figures_dir, "wgcna_module_trait.png")
    png(outfile, width = 10, height = max(6, nrow(mt_matrix) * 0.4), units = "in", res = 150)

    ht <- Heatmap(mt_matrix,
                  name = "Correlation",
                  col = col_fun,
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 9),
                  column_title = "WGCNA Module-Trait Correlations",
                  heatmap_legend_param = list(title = "r"))
    draw(ht)
    dev.off()
  } else {
    # Fallback: ggplot tile heatmap
    mt_long <- melt(as.data.table(mt_data), id.vars = names(mt_data)[1],
                    variable.name = "trait", value.name = "correlation")
    names(mt_long)[1] <- "module"

    p <- ggplot(mt_long, aes(x = trait, y = module, fill = correlation)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                          midpoint = 0, limits = c(-1, 1)) +
      labs(
        title = "WGCNA Module-Trait Correlations",
        x = "Trait",
        y = "Module",
        fill = "Correlation"
      ) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    outfile <- file.path(figures_dir, "wgcna_module_trait.png")
    ggsave(outfile, p, width = 10, height = max(6, nrow(mt_data) * 0.4), dpi = 150)
  }

  register_figure("wgcna_module_trait.png", "WGCNA Module-Trait Correlations",
                  "WGCNA", "Heatmap of module-trait correlation values")
  cat("  Generated: wgcna_module_trait.png\n")
}

create_wgcna_figures()

# ============================================================================
# 5. OVERVIEW - SAMPLE-LEVEL DASHBOARD
# ============================================================================

cat("=== Generating Overview Dashboard ===\n")

create_overview_dashboard <- function() {
  # Gather sample-level information from multiple sources
  dashboard_data <- data.frame(sample = character(), stringsAsFactors = FALSE)

  # Try to load QC summary
  qc_dir <- file.path(opt$results_dir, "02_qc")
  if (!dir.exists(qc_dir)) qc_dir <- file.path(opt$results_dir, "qc")

  multiqc_file <- list.files(qc_dir, pattern = "multiqc_general_stats",
                            full.names = TRUE, recursive = TRUE)[1]
  if (!is.na(multiqc_file) && file.exists(multiqc_file)) {
    qc_data <- safe_fread(multiqc_file)
    if (!is.null(qc_data)) {
      sample_col <- names(qc_data)[1]
      dashboard_data <- data.frame(sample = qc_data[[sample_col]], stringsAsFactors = FALSE)
      # Add QC pass/fail based on available metrics
      total_reads_col <- grep("total_reads|total_sequences|M Seqs", names(qc_data),
                             value = TRUE, ignore.case = TRUE)[1]
      if (!is.na(total_reads_col)) {
        dashboard_data$total_reads <- qc_data[[total_reads_col]]
        dashboard_data$qc_status <- ifelse(dashboard_data$total_reads > 1e6, "PASS", "FAIL")
      }
    }
  }

  # If no samples gathered yet, try from normalized counts
  if (nrow(dashboard_data) == 0) {
    norm_file <- file.path(opt$results_dir, "03_counting", "normalized", "vst_counts.tsv")
    if (!file.exists(norm_file)) {
      norm_file <- file.path(opt$results_dir, "counting", "normalized", "vst_counts.tsv")
    }
    if (file.exists(norm_file)) {
      counts_header <- names(fread(norm_file, nrows = 0))
      samples <- counts_header[-1]
      dashboard_data <- data.frame(sample = samples, stringsAsFactors = FALSE)
      dashboard_data$qc_status <- "PASS"  # default
    }
  }

  if (nrow(dashboard_data) == 0) {
    cat("  Cannot build overview dashboard - no sample data found.\n")
    return(NULL)
  }

  # Try to add ancestry information
  anc_dir <- file.path(opt$results_dir, "04_ancestry")
  if (!dir.exists(anc_dir)) anc_dir <- file.path(opt$results_dir, "ancestry")
  anc_file <- list.files(anc_dir, pattern = "proportions|ancestry_results",
                        full.names = TRUE)[1]
  if (!is.na(anc_file) && file.exists(anc_file)) {
    anc_data <- safe_fread(anc_file)
    if (!is.null(anc_data)) {
      pop_col <- grep("population|ancestry|category|predicted", names(anc_data),
                     value = TRUE, ignore.case = TRUE)[1]
      if (!is.na(pop_col)) {
        anc_merge <- data.frame(
          sample = anc_data[[1]],
          ancestry = anc_data[[pop_col]],
          stringsAsFactors = FALSE
        )
        dashboard_data <- merge(dashboard_data, anc_merge, by = "sample", all.x = TRUE)
      }
    }
  }

  # Try to add subtyping information
  sub_dir <- file.path(opt$results_dir, "18_molecular_subtyping")
  if (!dir.exists(sub_dir)) sub_dir <- file.path(opt$results_dir, "molecular_subtyping")
  sub_file <- list.files(sub_dir, pattern = "subtype|classification",
                        full.names = TRUE, recursive = TRUE)[1]
  if (!is.na(sub_file) && file.exists(sub_file)) {
    sub_data <- safe_fread(sub_file)
    if (!is.null(sub_data)) {
      sub_col <- grep("subtype|class|cluster", names(sub_data),
                     value = TRUE, ignore.case = TRUE)[1]
      if (!is.na(sub_col)) {
        sub_merge <- data.frame(
          sample = sub_data[[1]],
          subtype = sub_data[[sub_col]],
          stringsAsFactors = FALSE
        )
        dashboard_data <- merge(dashboard_data, sub_merge, by = "sample", all.x = TRUE)
      }
    }
  }

  # Try to add TMB information
  var_dir <- file.path(opt$results_dir, "05_variant_calling")
  if (!dir.exists(var_dir)) var_dir <- file.path(opt$results_dir, "variant_calling")
  tmb_file <- list.files(var_dir, pattern = "tmb|mutation_burden",
                        full.names = TRUE, recursive = TRUE)[1]
  if (!is.na(tmb_file) && file.exists(tmb_file)) {
    tmb_data <- safe_fread(tmb_file)
    if (!is.null(tmb_data)) {
      tmb_col <- grep("tmb|burden|mutations_per_mb", names(tmb_data),
                     value = TRUE, ignore.case = TRUE)[1]
      if (!is.na(tmb_col)) {
        tmb_merge <- data.frame(
          sample = tmb_data[[1]],
          tmb = tmb_data[[tmb_col]],
          stringsAsFactors = FALSE
        )
        dashboard_data <- merge(dashboard_data, tmb_merge, by = "sample", all.x = TRUE)
      }
    }
  }

  # Build dashboard figure using patchwork
  plots <- list()

  # QC status bar
  if ("qc_status" %in% names(dashboard_data)) {
    qc_summary <- as.data.frame(table(dashboard_data$qc_status))
    names(qc_summary) <- c("Status", "Count")
    p_qc <- ggplot(qc_summary, aes(x = Status, y = Count, fill = Status)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("PASS" = "#4DAF4A", "FAIL" = "#E41A1C")) +
      labs(title = "QC Status", x = "", y = "Samples") +
      theme_minimal() +
      theme(legend.position = "none")
    plots[["qc"]] <- p_qc
  }

  # Ancestry bar
  if ("ancestry" %in% names(dashboard_data)) {
    anc_summary <- as.data.frame(table(dashboard_data$ancestry))
    names(anc_summary) <- c("Ancestry", "Count")
    p_anc <- ggplot(anc_summary, aes(x = reorder(Ancestry, -Count), y = Count, fill = Ancestry)) +
      geom_bar(stat = "identity") +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "Ancestry Distribution", x = "", y = "Samples") +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
    plots[["ancestry"]] <- p_anc
  }

  # Subtype bar
  if ("subtype" %in% names(dashboard_data)) {
    sub_summary <- as.data.frame(table(dashboard_data$subtype))
    names(sub_summary) <- c("Subtype", "Count")
    p_sub <- ggplot(sub_summary, aes(x = reorder(Subtype, -Count), y = Count, fill = Subtype)) +
      geom_bar(stat = "identity") +
      scale_fill_brewer(palette = "Set3") +
      labs(title = "Molecular Subtypes", x = "", y = "Samples") +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
    plots[["subtype"]] <- p_sub
  }

  # TMB distribution
  if ("tmb" %in% names(dashboard_data) && is.numeric(dashboard_data$tmb)) {
    p_tmb <- ggplot(dashboard_data, aes(x = tmb)) +
      geom_histogram(bins = 30, fill = "#984EA3", alpha = 0.7) +
      geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
      labs(title = "Tumor Mutation Burden", x = "TMB (mutations/Mb)", y = "Samples") +
      theme_minimal()
    plots[["tmb"]] <- p_tmb
  }

  if (length(plots) == 0) {
    cat("  Not enough data for overview dashboard.\n")
    return(NULL)
  }

  # Combine with patchwork
  combined_plot <- wrap_plots(plots, ncol = min(2, length(plots))) +
    plot_annotation(
      title = "Sample Overview Dashboard",
      subtitle = sprintf("N = %d samples", nrow(dashboard_data)),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    )

  outfile <- file.path(figures_dir, "overview_dashboard.png")
  n_rows <- ceiling(length(plots) / 2)
  ggsave(outfile, combined_plot, width = 12, height = 5 * n_rows, dpi = 150)
  register_figure("overview_dashboard.png", "Sample Overview Dashboard",
                  "Overview", "Combined dashboard showing QC, ancestry, subtype, and TMB")
  cat("  Generated: overview_dashboard.png\n")

  # Save dashboard data
  fwrite(dashboard_data, file.path(opt$output_dir, "dashboard_data.tsv"), sep = "\t")
}

create_overview_dashboard()

# ============================================================================
# CREATE HTML FIGURE INDEX
# ============================================================================

cat("=== Creating HTML Figure Index ===\n")

create_html_index <- function() {
  if (nrow(figure_registry) == 0) {
    cat("  No figures generated, skipping HTML index.\n")
    return(NULL)
  }

  # Group by category
  categories <- unique(figure_registry$category)

  # Build HTML
  html_content <- paste0(
    '<!DOCTYPE html>\n',
    '<html lang="en">\n',
    '<head>\n',
    '  <meta charset="UTF-8">\n',
    '  <meta name="viewport" content="width=device-width, initial-scale=1.0">\n',
    '  <title>Pipeline Visualization Index</title>\n',
    '  <style>\n',
    '    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;\n',
    '           max-width: 1200px; margin: 0 auto; padding: 20px; background: #f5f5f5; }\n',
    '    h1 { color: #333; border-bottom: 2px solid #4CAF50; padding-bottom: 10px; }\n',
    '    h2 { color: #555; margin-top: 30px; }\n',
    '    .figure-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));\n',
    '                   gap: 20px; margin-top: 15px; }\n',
    '    .figure-card { background: white; border-radius: 8px; padding: 15px;\n',
    '                   box-shadow: 0 2px 4px rgba(0,0,0,0.1); transition: transform 0.2s; }\n',
    '    .figure-card:hover { transform: translateY(-2px); box-shadow: 0 4px 8px rgba(0,0,0,0.15); }\n',
    '    .figure-card img { width: 100%; height: auto; border-radius: 4px; border: 1px solid #eee; }\n',
    '    .figure-card h3 { margin: 10px 0 5px; font-size: 14px; color: #333; }\n',
    '    .figure-card p { margin: 0; font-size: 12px; color: #666; }\n',
    '    .timestamp { color: #999; font-size: 12px; margin-top: 30px; }\n',
    '  </style>\n',
    '</head>\n',
    '<body>\n',
    '  <h1>Cancer RNA-seq Pipeline - Figure Index</h1>\n',
    sprintf('  <p>Generated: %s | Total figures: %d</p>\n',
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), nrow(figure_registry))
  )

  for (cat_name in categories) {
    cat_figures <- figure_registry[figure_registry$category == cat_name, ]
    html_content <- paste0(html_content, sprintf('  <h2>%s</h2>\n  <div class="figure-grid">\n', cat_name))

    for (i in seq_len(nrow(cat_figures))) {
      fig <- cat_figures[i, ]
      fig_path <- file.path("../figures", fig$filename)
      html_content <- paste0(html_content, sprintf(
        '    <div class="figure-card">\n      <a href="%s"><img src="%s" alt="%s"></a>\n      <h3>%s</h3>\n      <p>%s</p>\n    </div>\n',
        fig_path, fig_path, fig$title, fig$title, fig$description
      ))
    }

    html_content <- paste0(html_content, '  </div>\n')
  }

  html_content <- paste0(html_content,
    sprintf('  <p class="timestamp">Report generated by cancer RNA-seq pipeline visualization step.</p>\n'),
    '</body>\n</html>\n'
  )

  outfile <- file.path(html_dir, "figure_index.html")
  writeLines(html_content, outfile)
  cat(sprintf("  HTML index: %s\n", outfile))
}

create_html_index()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat(sprintf("\n=== Visualization Complete ===\n"))
cat(sprintf("  Figures generated: %d\n", nrow(figure_registry)))
cat(sprintf("  Output directory: %s\n", opt$output_dir))
cat(sprintf("  HTML index: %s\n", file.path(html_dir, "figure_index.html")))
