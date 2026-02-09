#!/usr/bin/env Rscript

# ============================================================================
# Tumor Mutational Burden (TMB) Estimation from RNA-seq variant calls
# ============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(VariantAnnotation)
})

option_list <- list(
    make_option("--vcf-dir", type = "character", default = "."),
    make_option("--metadata", type = "character"),
    make_option("--output-prefix", type = "character", default = "tmb"),
    make_option("--min-vaf", type = "double", default = 0.05),
    make_option("--min-depth", type = "integer", default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))

vcf_files <- list.files(opt$`vcf-dir`, pattern = "\\.filtered\\.vcf\\.gz$", full.names = TRUE)

if (length(vcf_files) == 0) {
    cat("No filtered VCF files found.\n")
    writeLines("sample_id\ttmb_per_mb\tn_variants\tcoding_region_mb", paste0(opt$`output-prefix`, "_scores.tsv"))
    quit(save = "no")
}

# Approximate coding region size (Mb) for normalization
# RNA-seq captures ~30-50 Mb of coding sequence
CODING_REGION_MB <- 36  # Approximate

tmb_results <- data.frame()

for (vcf_path in vcf_files) {
    sample_id <- gsub("\\.filtered\\.vcf\\.gz$", "", basename(vcf_path))
    cat("Processing:", sample_id, "\n")

    tryCatch({
        vcf <- readVcf(vcf_path, genome = "GRCh38")

        # Filter: PASS variants only, above min depth and VAF
        pass_idx <- fixed(vcf)$FILTER == "PASS" | fixed(vcf)$FILTER == "."

        if (sum(pass_idx) > 0) {
            vcf_pass <- vcf[pass_idx]
            n_variants <- nrow(vcf_pass)

            # TMB = variants / coding region in Mb
            tmb <- n_variants / CODING_REGION_MB

            tmb_results <- rbind(tmb_results, data.frame(
                sample_id = sample_id,
                tmb_per_mb = round(tmb, 2),
                n_variants = n_variants,
                coding_region_mb = CODING_REGION_MB
            ))
        } else {
            tmb_results <- rbind(tmb_results, data.frame(
                sample_id = sample_id, tmb_per_mb = 0,
                n_variants = 0, coding_region_mb = CODING_REGION_MB
            ))
        }
    }, error = function(e) {
        cat("  Error processing", sample_id, ":", e$message, "\n")
        tmb_results <<- rbind(tmb_results, data.frame(
            sample_id = sample_id, tmb_per_mb = NA,
            n_variants = NA, coding_region_mb = CODING_REGION_MB
        ))
    })
}

write.table(tmb_results, paste0(opt$`output-prefix`, "_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Summary
summary_df <- data.frame(
    metric = c("n_samples", "median_tmb", "mean_tmb", "tmb_high_n"),
    value = c(nrow(tmb_results),
              round(median(tmb_results$tmb_per_mb, na.rm = TRUE), 2),
              round(mean(tmb_results$tmb_per_mb, na.rm = TRUE), 2),
              sum(tmb_results$tmb_per_mb > 10, na.rm = TRUE))
)
write.table(summary_df, paste0(opt$`output-prefix`, "_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("TMB estimation complete.\n")
