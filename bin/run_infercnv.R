#!/usr/bin/env Rscript

# ============================================================================
# InferCNV - Copy Number Inference from RNA-seq
# ============================================================================

suppressPackageStartupMessages({
    library(infercnv)
    library(optparse)
})

option_list <- list(
    make_option("--counts", type = "character"),
    make_option("--gene-order", type = "character"),
    make_option("--annotations", type = "character"),
    make_option("--output-dir", type = "character", default = "infercnv_output"),
    make_option("--num-threads", type = "integer", default = 4)
)
opt <- parse_args(OptionParser(option_list = option_list))

counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
gene_order <- read.delim(opt$`gene-order`, header = FALSE, row.names = 1)
annotations <- read.delim(opt$annotations, header = FALSE, row.names = 1)

# Create InferCNV object
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = as.matrix(round(counts)),
    annotations_file = opt$annotations,
    gene_order_file = opt$`gene-order`,
    ref_group_names = c("reference")
)

# Run InferCNV
infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = opt$`output-dir`,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    num_threads = opt$`num-threads`
)

cat("InferCNV analysis complete.\n")
