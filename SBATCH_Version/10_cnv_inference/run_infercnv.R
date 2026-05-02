#!/usr/bin/env Rscript
# ============================================================================
# InferCNV: Copy Number Variation Inference from RNA-Seq
# ============================================================================
# Infers large-scale chromosomal copy number variations by comparing expression
# patterns in tumor samples against a set of reference (normal) samples.
#
# Usage:
#   Rscript run_infercnv.R \
#     --counts_matrix raw_counts_matrix.tsv \
#     --gene_order_file gene_positions.bed \
#     --annotations_file annotations.txt \
#     --ref_group reference \
#     --output_dir /path/to/output \
#     --num_threads 16
# ============================================================================

suppressPackageStartupMessages({
    library(infercnv)
    library(optparse)
})

# --- Parse command-line arguments ---
option_list <- list(
    make_option(c("--counts_matrix"), type = "character", default = NULL,
                help = "Path to raw counts matrix (genes x samples, tab-delimited)"),
    make_option(c("--gene_order_file"), type = "character", default = NULL,
                help = "Path to gene order/position file (gene, chr, start, stop)"),
    make_option(c("--annotations_file"), type = "character", default = NULL,
                help = "Path to annotations file (sample_id, group)"),
    make_option(c("--ref_group"), type = "character", default = "reference",
                help = "Name of the reference group in annotations file [default: %default]"),
    make_option(c("--output_dir"), type = "character", default = "infercnv_output",
                help = "Output directory [default: %default]"),
    make_option(c("--num_threads"), type = "integer", default = 4,
                help = "Number of threads [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# --- Validate required arguments ---
if (is.null(opt$counts_matrix)) stop("--counts_matrix is required")
if (is.null(opt$gene_order_file)) stop("--gene_order_file is required")
if (is.null(opt$annotations_file)) stop("--annotations_file is required")

cat("=== InferCNV Analysis ===\n")
cat("Counts matrix:", opt$counts_matrix, "\n")
cat("Gene order file:", opt$gene_order_file, "\n")
cat("Annotations file:", opt$annotations_file, "\n")
cat("Reference group:", opt$ref_group, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Threads:", opt$num_threads, "\n")

# --- Create output directory ---
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Create InferCNV object ---
cat("\nCreating InferCNV object...\n")

infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = opt$counts_matrix,
    annotations_file = opt$annotations_file,
    gene_order_file = opt$gene_order_file,
    ref_group_names = c(opt$ref_group)
)

# --- Run InferCNV ---
cat("Running InferCNV with HMM...\n")

infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = opt$output_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    num_threads = opt$num_threads,
    analysis_mode = "subclusters",
    tumor_subcluster_partition_method = "random_trees",
    output_format = "pdf",
    plot_steps = FALSE
)

# --- Save results ---
cat("\nSaving InferCNV results...\n")

# Save the final object
saveRDS(infercnv_obj, file = file.path(opt$output_dir, "infercnv_obj.rds"))

# Export HMM predictions if available
hmm_results_file <- file.path(opt$output_dir, "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat")
if (file.exists(hmm_results_file)) {
    cat("HMM predictions saved to:", hmm_results_file, "\n")
}

# Export observations matrix
observations_file <- file.path(opt$output_dir, "infercnv.observations.txt")
if (file.exists(observations_file)) {
    cat("Observations matrix saved to:", observations_file, "\n")
}

cat("\n=== InferCNV analysis complete ===\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Key output files:\n")
cat("  - infercnv.pdf (heatmap)\n")
cat("  - infercnv_obj.rds (R object)\n")
cat("  - HMM CNV predictions\n")
cat("  - infercnv.observations.txt\n")
