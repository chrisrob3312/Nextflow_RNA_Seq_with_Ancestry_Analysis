#!/usr/bin/env Rscript

# ============================================================================
# GSVA - Gene Set Variation Analysis (sample-level pathway scores)
# ============================================================================

suppressPackageStartupMessages({
    library(GSVA)
    library(msigdbr)
    library(optparse)
    library(ggplot2)
    library(pheatmap)
})

option_list <- list(
    make_option("--expression", type = "character"),
    make_option("--metadata", type = "character"),
    make_option("--databases", type = "character", default = "HALLMARK,IMMUNESIGDB"),
    make_option("--species", type = "character", default = "Homo sapiens"),
    make_option("--min-size", type = "integer", default = 15),
    make_option("--max-size", type = "integer", default = 500),
    make_option("--output-scores", type = "character", default = "gsva_scores.tsv"),
    make_option("--output-dir", type = "character", default = "gsva_results"),
    make_option("--plot-dir", type = "character", default = "gsva_plots")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)

expr <- as.matrix(read.delim(opt$expression, row.names = 1, check.names = FALSE))
metadata <- read.delim(opt$metadata, check.names = FALSE)
rownames(metadata) <- metadata$sample_id

db_names <- trimws(unlist(strsplit(opt$databases, ",")))
all_scores <- list()

for (db_name in db_names) {
    cat("Running GSVA for:", db_name, "\n")

    msig <- tryCatch({
        if (db_name == "HALLMARK") msigdbr(species = opt$species, category = "H")
        else if (db_name == "IMMUNESIGDB") msigdbr(species = opt$species, category = "C7", subcategory = "IMMUNESIGDB")
        else if (db_name == "KEGG") msigdbr(species = opt$species, category = "C2", subcategory = "CP:KEGG")
        else if (db_name == "REACTOME") msigdbr(species = opt$species, category = "C2", subcategory = "CP:REACTOME")
        else if (db_name == "GO_BP") msigdbr(species = opt$species, category = "C5", subcategory = "GO:BP")
        else NULL
    }, error = function(e) NULL)

    if (is.null(msig)) next

    gene_sets <- split(msig$gene_symbol, msig$gs_name)
    gene_sets <- gene_sets[sapply(gene_sets, length) >= opt$`min-size` &
                           sapply(gene_sets, length) <= opt$`max-size`]

    if (length(gene_sets) == 0) next

    gsva_res <- gsva(expr, gene_sets, method = "gsva", kcdf = "Gaussian",
                     min.sz = opt$`min-size`, max.sz = opt$`max-size`, verbose = FALSE)

    write.table(gsva_res, file.path(opt$`output-dir`, paste0("gsva_", db_name, ".tsv")),
                sep = "\t", quote = FALSE)
    all_scores[[db_name]] <- gsva_res

    # Heatmap of top variable pathways
    tryCatch({
        vars <- apply(gsva_res, 1, var)
        top_pathways <- names(head(sort(vars, decreasing = TRUE), 30))
        annot <- metadata[colnames(gsva_res), intersect(c("cytomolecular_subgroup", "relapse_status"), colnames(metadata)), drop = FALSE]

        pdf(file.path(opt$`plot-dir`, paste0("gsva_", db_name, "_heatmap.pdf")), width = 14, height = 10)
        pheatmap(gsva_res[top_pathways, ], scale = "none",
                 annotation_col = annot,
                 main = paste("GSVA Scores -", db_name),
                 fontsize_row = 7, fontsize_col = 8)
        dev.off()
    }, error = function(e) cat("Heatmap error:", e$message, "\n"))
}

# Combined scores output
if (length(all_scores) > 0) {
    combined <- do.call(rbind, all_scores)
    write.table(combined, opt$`output-scores`, sep = "\t", quote = FALSE)
}

cat("GSVA analysis complete.\n")
