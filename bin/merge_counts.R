#!/usr/bin/env Rscript

# ============================================================================
# Merge individual featureCounts outputs into a count matrix
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
})

option_list <- list(
    make_option("--input-dir", type = "character", default = "."),
    make_option("--pattern", type = "character", default = "*.featureCounts.txt"),
    make_option("--output-matrix", type = "character", default = "raw_count_matrix.tsv"),
    make_option("--output-genes", type = "character", default = "gene_annotations.tsv"),
    make_option("--output-stats", type = "character", default = "count_summary_stats.tsv"),
    make_option("--min-counts", type = "integer", default = 10),
    make_option("--min-samples", type = "integer", default = 3)
)
opt <- parse_args(OptionParser(option_list = option_list))

# Find featureCounts files
files <- list.files(opt$`input-dir`, pattern = "\\.featureCounts\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("No featureCounts files found")

cat("Found", length(files), "featureCounts files\n")

# Read first file to get gene info
first <- read.delim(files[1], comment.char = "#", check.names = FALSE)
gene_info <- first[, 1:6]
colnames(gene_info) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

# Build count matrix
count_matrix <- matrix(0, nrow = nrow(first), ncol = length(files))
rownames(count_matrix) <- first$Geneid

sample_names <- c()
for (i in seq_along(files)) {
    dat <- read.delim(files[i], comment.char = "#", check.names = FALSE)
    # The count column is the last column (7th)
    count_matrix[, i] <- dat[, ncol(dat)]
    # Extract sample name from column header or filename
    sname <- colnames(dat)[ncol(dat)]
    sname <- gsub("\\.sorted\\.markdup\\.bam$", "", basename(sname))
    sname <- gsub("\\.bam$", "", sname)
    sample_names <- c(sample_names, sname)
}
colnames(count_matrix) <- sample_names

# Summary stats
stats <- data.frame(
    sample_id = sample_names,
    total_counts = colSums(count_matrix),
    genes_detected = colSums(count_matrix > 0),
    median_counts = apply(count_matrix, 2, median),
    stringsAsFactors = FALSE
)

# Write outputs
write.table(count_matrix, opt$`output-matrix`, sep = "\t", quote = FALSE)
write.table(gene_info, opt$`output-genes`, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(stats, opt$`output-stats`, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Count matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
cat("Written to:", opt$`output-matrix`, "\n")
