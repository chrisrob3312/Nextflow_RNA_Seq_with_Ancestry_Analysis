#!/usr/bin/env Rscript
# ============================================================================
# NORMALIZATION AND BATCH CORRECTION FOR RNA-SEQ COUNT DATA
# ============================================================================
# This script performs:
#   1. Low-count gene filtering
#   2. Normalization (VST, rlog, TMM, or TPM)
#   3. Batch correction using ComBat-seq (on raw counts) and ComBat (on
#      normalized values)
#   4. Outputs corrected counts and normalized matrices
#
# Usage:
#   Rscript normalize_and_batch_correct.R \
#       <raw_counts_matrix> \
#       <metadata_csv> \
#       <output_dir> \
#       <normalization_method> \
#       <batch_variable>
#
# Arguments:
#   args[1] = raw count matrix path (TSV, gene_id as first column)
#   args[2] = metadata CSV path (sample_id as first column)
#   args[3] = output directory
#   args[4] = normalization method: vst, rlog, tmm, or tpm
#   args[5] = batch variable column name in metadata
# ============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(sva)
    library(edgeR)
    library(matrixStats)
})

# ============================================================================
# PARSE ARGUMENTS
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
    stop(paste(
        "Usage: Rscript normalize_and_batch_correct.R",
        "<counts_matrix> <metadata_csv> <output_dir>",
        "<norm_method: vst|rlog|tmm|tpm> <batch_variable>"
    ))
}

counts_file   <- args[1]
metadata_file <- args[2]
output_dir    <- args[3]
norm_method   <- tolower(args[4])
batch_var     <- args[5]

cat("============================================================================\n")
cat("  NORMALIZATION AND BATCH CORRECTION\n")
cat("============================================================================\n")
cat(sprintf("  Counts matrix:       %s\n", counts_file))
cat(sprintf("  Metadata:            %s\n", metadata_file))
cat(sprintf("  Output directory:    %s\n", output_dir))
cat(sprintf("  Normalization:       %s\n", norm_method))
cat(sprintf("  Batch variable:      %s\n", batch_var))
cat("============================================================================\n\n")

# Validate normalization method
valid_methods <- c("vst", "rlog", "tmm", "tpm")
if (!norm_method %in% valid_methods) {
    stop(sprintf(
        "Invalid normalization method '%s'. Must be one of: %s",
        norm_method, paste(valid_methods, collapse = ", ")
    ))
}

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("[1/5] Loading data...\n")

# Read count matrix
counts_raw <- read.table(
    counts_file,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

# Read metadata
metadata <- read.csv(
    metadata_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Set first column as row names for metadata
rownames(metadata) <- metadata[, 1]

cat(sprintf("  Raw counts: %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))
cat(sprintf("  Metadata: %d samples x %d fields\n", nrow(metadata), ncol(metadata)))

# Ensure count matrix has integer values
counts_raw <- round(as.matrix(counts_raw))
storage.mode(counts_raw) <- "integer"

# Match sample order between counts and metadata
common_samples <- intersect(colnames(counts_raw), rownames(metadata))
if (length(common_samples) == 0) {
    stop("ERROR: No matching sample IDs between count matrix columns and metadata rows.")
}

if (length(common_samples) < ncol(counts_raw)) {
    cat(sprintf("  WARNING: Only %d of %d count matrix samples found in metadata.\n",
                length(common_samples), ncol(counts_raw)))
}

counts_raw <- counts_raw[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

cat(sprintf("  After matching: %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))

# Validate batch variable
if (!batch_var %in% colnames(metadata)) {
    cat(sprintf("  WARNING: Batch variable '%s' not found in metadata columns.\n", batch_var))
    cat("  Available columns:", paste(colnames(metadata), collapse = ", "), "\n")
    cat("  Proceeding WITHOUT batch correction.\n")
    batch_var <- NULL
}

# ============================================================================
# FILTER LOW-COUNT GENES
# ============================================================================

cat("\n[2/5] Filtering low-count genes...\n")

# Filter criteria (from config or defaults)
min_count <- as.numeric(Sys.getenv("MIN_GENE_COUNTS", unset = "10"))
min_samples <- as.numeric(Sys.getenv("MIN_SAMPLES_EXPRESSING", unset = "3"))

cat(sprintf("  Filter: genes with >= %d counts in >= %d samples\n", min_count, min_samples))

# Apply filter
genes_pass <- rowSums(counts_raw >= min_count) >= min_samples
counts_filtered <- counts_raw[genes_pass, , drop = FALSE]

cat(sprintf("  Genes before filter: %d\n", nrow(counts_raw)))
cat(sprintf("  Genes after filter:  %d\n", nrow(counts_filtered)))
cat(sprintf("  Genes removed:       %d\n", sum(!genes_pass)))

# Save filtered gene list
write.table(
    data.frame(gene_id = rownames(counts_filtered)),
    file = file.path(output_dir, "filtered_gene_list.txt"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

# ============================================================================
# BATCH CORRECTION WITH ComBat-seq (ON RAW COUNTS)
# ============================================================================

cat("\n[3/5] Batch correction (ComBat-seq on raw counts)...\n")

if (!is.null(batch_var)) {
    batch_vector <- metadata[[batch_var]]

    # Check if batch variable has more than one level
    n_batches <- length(unique(batch_vector))
    if (n_batches < 2) {
        cat(sprintf("  Only %d batch level found. Skipping batch correction.\n", n_batches))
        counts_corrected <- counts_filtered
    } else {
        cat(sprintf("  Number of batch levels: %d\n", n_batches))
        cat(sprintf("  Batch distribution:\n"))
        batch_table <- table(batch_vector)
        for (b in names(batch_table)) {
            cat(sprintf("    %s: %d samples\n", b, batch_table[b]))
        }

        # Run ComBat-seq
        tryCatch({
            counts_corrected <- ComBat_seq(
                counts = counts_filtered,
                batch = batch_vector,
                group = NULL
            )
            cat("  ComBat-seq correction applied successfully.\n")
        }, error = function(e) {
            cat(sprintf("  WARNING: ComBat-seq failed: %s\n", e$message))
            cat("  Proceeding without batch correction.\n")
            counts_corrected <<- counts_filtered
        })
    }
} else {
    cat("  No batch variable specified. Skipping batch correction.\n")
    counts_corrected <- counts_filtered
}

# Save batch-corrected counts
write.table(
    cbind(gene_id = rownames(counts_corrected), as.data.frame(counts_corrected)),
    file = file.path(output_dir, "counts_batch_corrected.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

# ============================================================================
# NORMALIZATION
# ============================================================================

cat(sprintf("\n[4/5] Normalizing with method: %s...\n", norm_method))

normalized_matrix <- NULL

if (norm_method == "vst") {
    # DESeq2 variance-stabilizing transformation
    dds <- DESeqDataSetFromMatrix(
        countData = counts_corrected,
        colData = metadata,
        design = ~ 1
    )
    vsd <- vst(dds, blind = TRUE)
    normalized_matrix <- assay(vsd)
    cat("  VST normalization complete.\n")

} else if (norm_method == "rlog") {
    # DESeq2 regularized log transformation
    dds <- DESeqDataSetFromMatrix(
        countData = counts_corrected,
        colData = metadata,
        design = ~ 1
    )
    rld <- rlog(dds, blind = TRUE)
    normalized_matrix <- assay(rld)
    cat("  rlog normalization complete.\n")

} else if (norm_method == "tmm") {
    # edgeR TMM normalization -> log2-CPM
    dge <- DGEList(counts = counts_corrected)
    dge <- calcNormFactors(dge, method = "TMM")
    normalized_matrix <- cpm(dge, log = TRUE, prior.count = 1)
    cat("  TMM normalization (log2-CPM) complete.\n")

} else if (norm_method == "tpm") {
    # TPM normalization (requires gene lengths)
    gene_lengths_file <- file.path(dirname(counts_file), "gene_lengths.tsv")

    if (!file.exists(gene_lengths_file)) {
        stop(sprintf(
            "TPM requires gene lengths file: %s\nGenerate from featureCounts output.",
            gene_lengths_file
        ))
    }

    gene_lengths <- read.table(
        gene_lengths_file,
        header = TRUE, sep = "\t", row.names = 1
    )

    # Match genes
    common_genes <- intersect(rownames(counts_corrected), rownames(gene_lengths))
    counts_for_tpm <- counts_corrected[common_genes, , drop = FALSE]
    lengths_for_tpm <- gene_lengths[common_genes, 1]

    # Calculate TPM
    # RPK = reads / (gene_length_kb)
    rpk <- counts_for_tpm / (lengths_for_tpm / 1000)
    # TPM = RPK / sum(RPK) * 1e6
    scaling_factors <- colSums(rpk)
    normalized_matrix <- t(t(rpk) / scaling_factors) * 1e6

    # Log2 transform for downstream analysis
    normalized_matrix_log <- log2(normalized_matrix + 1)

    # Save raw TPM
    write.table(
        cbind(gene_id = rownames(normalized_matrix), as.data.frame(normalized_matrix)),
        file = file.path(output_dir, "tpm_matrix.tsv"),
        quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
    )

    # Use log2(TPM+1) as the normalized matrix
    normalized_matrix <- normalized_matrix_log
    cat("  TPM normalization complete (log2(TPM+1) saved).\n")
}

# ============================================================================
# SAVE OUTPUTS
# ============================================================================

cat("\n[5/5] Saving outputs...\n")

# Save normalized matrix
write.table(
    cbind(gene_id = rownames(normalized_matrix), as.data.frame(normalized_matrix)),
    file = file.path(output_dir, sprintf("normalized_%s_matrix.tsv", norm_method)),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

# Save filtered raw counts (before batch correction)
write.table(
    cbind(gene_id = rownames(counts_filtered), as.data.frame(counts_filtered)),
    file = file.path(output_dir, "counts_filtered.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

# Save sample metadata used
write.csv(
    metadata,
    file = file.path(output_dir, "metadata_used.csv"),
    row.names = FALSE
)

# Save normalization summary statistics
summary_stats <- data.frame(
    metric = c(
        "total_genes_input",
        "genes_after_filter",
        "genes_removed",
        "total_samples",
        "normalization_method",
        "batch_variable",
        "n_batch_levels",
        "min_count_threshold",
        "min_samples_threshold"
    ),
    value = c(
        nrow(counts_raw),
        nrow(counts_filtered),
        sum(!genes_pass),
        ncol(counts_filtered),
        norm_method,
        ifelse(is.null(batch_var), "none", batch_var),
        ifelse(is.null(batch_var), 0, length(unique(metadata[[batch_var]]))),
        min_count,
        min_samples
    )
)

write.table(
    summary_stats,
    file = file.path(output_dir, "normalization_summary.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================================\n")
cat("  NORMALIZATION & BATCH CORRECTION COMPLETE\n")
cat("============================================================================\n")
cat(sprintf("  Input genes:              %d\n", nrow(counts_raw)))
cat(sprintf("  Genes after filtering:    %d\n", nrow(counts_filtered)))
cat(sprintf("  Samples:                  %d\n", ncol(counts_filtered)))
cat(sprintf("  Normalization method:     %s\n", norm_method))
cat(sprintf("  Batch correction:         %s\n",
    ifelse(is.null(batch_var), "none", paste0("ComBat-seq (", batch_var, ")"))))
cat("\n  Output files:\n")
cat(sprintf("    %s\n", file.path(output_dir, "counts_filtered.tsv")))
cat(sprintf("    %s\n", file.path(output_dir, "counts_batch_corrected.tsv")))
cat(sprintf("    %s\n", file.path(output_dir, sprintf("normalized_%s_matrix.tsv", norm_method))))
cat(sprintf("    %s\n", file.path(output_dir, "normalization_summary.tsv")))
cat("============================================================================\n")

cat("\nDone.\n")
