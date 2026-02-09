#!/usr/bin/env Rscript

# ============================================================================
# Pathway Enrichment - ORA (clusterProfiler) + GSEA (fgsea) + MSigDB
# Databases: GO_BP, GO_MF, KEGG, Reactome, Hallmark, ImmuneSigDB
# ============================================================================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(fgsea)
    library(msigdbr)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(optparse)
    library(ggplot2)
})

option_list <- list(
    make_option("--de-results-dir", type = "character"),
    make_option("--counts", type = "character"),
    make_option("--databases", type = "character", default = "GO_BP,GO_MF,KEGG,REACTOME,HALLMARK,IMMUNESIGDB"),
    make_option("--species", type = "character", default = "Homo sapiens"),
    make_option("--min-size", type = "integer", default = 15),
    make_option("--max-size", type = "integer", default = 500),
    make_option("--padj-threshold", type = "double", default = 0.05),
    make_option("--output-dir", type = "character", default = "pathway_results"),
    make_option("--plot-dir", type = "character", default = "pathway_plots")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(file.path(opt$`output-dir`, "ora"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$`output-dir`, "gsea"), recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- Load gene sets from MSigDB ----
db_names <- trimws(unlist(strsplit(opt$databases, ",")))

load_genesets <- function(db_name, species) {
    if (db_name == "GO_BP") {
        msigdbr(species = species, category = "C5", subcategory = "GO:BP")
    } else if (db_name == "GO_MF") {
        msigdbr(species = species, category = "C5", subcategory = "GO:MF")
    } else if (db_name == "KEGG") {
        msigdbr(species = species, category = "C2", subcategory = "CP:KEGG")
    } else if (db_name == "REACTOME") {
        msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME")
    } else if (db_name == "HALLMARK") {
        msigdbr(species = species, category = "H")
    } else if (db_name == "IMMUNESIGDB") {
        msigdbr(species = species, category = "C7", subcategory = "IMMUNESIGDB")
    } else {
        NULL
    }
}

# ---- Find DE result files ----
de_files <- list.files(opt$`de-results-dir`, pattern = "\\.tsv$", full.names = TRUE)
de_files <- de_files[!grepl("summary", de_files)]

cat("Found", length(de_files), "DE result files\n")

# ---- Run enrichment for each DE result ----
for (de_file in de_files) {
    contrast_name <- gsub("\\.tsv$", "", basename(de_file))
    cat("\n==== Enrichment for:", contrast_name, "====\n")

    de_res <- read.delim(de_file, check.names = FALSE)
    if (!"log2FoldChange" %in% colnames(de_res) && "logFC" %in% colnames(de_res)) {
        de_res$log2FoldChange <- de_res$logFC
        de_res$padj <- de_res$adj.P.Val
    }

    # Ranked gene list for GSEA (by signed -log10 p-value)
    de_res$rank_score <- -log10(pmax(de_res$pvalue, 1e-300)) * sign(de_res$log2FoldChange)
    de_res <- de_res[!is.na(de_res$rank_score), ]
    gene_list <- setNames(de_res$rank_score, rownames(de_res))
    gene_list <- sort(gene_list, decreasing = TRUE)

    # Significant genes for ORA
    sig_up <- rownames(de_res)[de_res$padj < opt$`padj-threshold` & de_res$log2FoldChange > 0]
    sig_down <- rownames(de_res)[de_res$padj < opt$`padj-threshold` & de_res$log2FoldChange < 0]
    background <- rownames(de_res)

    for (db_name in db_names) {
        cat("  Database:", db_name, "\n")
        msig <- load_genesets(db_name, opt$species)
        if (is.null(msig) || nrow(msig) == 0) next

        # Convert to list for fgsea
        pathways_list <- split(msig$gene_symbol, msig$gs_name)

        # GSEA with fgsea
        tryCatch({
            gsea_res <- fgsea(
                pathways = pathways_list,
                stats = gene_list,
                minSize = opt$`min-size`,
                maxSize = opt$`max-size`,
                nproc = 1
            )
            gsea_res <- gsea_res[order(gsea_res$padj), ]
            write.table(as.data.frame(gsea_res[, c("pathway", "pval", "padj", "ES", "NES", "size")]),
                        file.path(opt$`output-dir`, "gsea", paste0(contrast_name, "_", db_name, "_gsea.tsv")),
                        sep = "\t", quote = FALSE, row.names = FALSE)

            # Plot top pathways
            if (sum(gsea_res$padj < opt$`padj-threshold`) > 0) {
                top_paths <- head(gsea_res[gsea_res$padj < opt$`padj-threshold`, ], 20)
                p <- ggplot(top_paths, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
                    geom_col() +
                    coord_flip() +
                    scale_fill_gradient(low = "red", high = "blue") +
                    labs(title = paste(contrast_name, "-", db_name, "GSEA"),
                         x = "", y = "Normalized Enrichment Score") +
                    theme_bw() +
                    theme(axis.text.y = element_text(size = 7))
                ggsave(file.path(opt$`plot-dir`, paste0(contrast_name, "_", db_name, "_gsea_barplot.pdf")),
                       p, width = 12, height = 8)
            }
        }, error = function(e) cat("    GSEA error:", e$message, "\n"))

        # ORA for up-regulated genes
        if (length(sig_up) > 5) {
            tryCatch({
                term2gene <- msig[, c("gs_name", "gene_symbol")]
                ora_up <- enricher(gene = sig_up, universe = background,
                                    TERM2GENE = term2gene,
                                    pAdjustMethod = "BH", pvalueCutoff = opt$`padj-threshold`,
                                    minGSSize = opt$`min-size`, maxGSSize = opt$`max-size`)
                if (!is.null(ora_up) && nrow(as.data.frame(ora_up)) > 0) {
                    write.table(as.data.frame(ora_up),
                                file.path(opt$`output-dir`, "ora", paste0(contrast_name, "_", db_name, "_ora_up.tsv")),
                                sep = "\t", quote = FALSE, row.names = FALSE)
                }
            }, error = function(e) cat("    ORA (up) error:", e$message, "\n"))
        }

        # ORA for down-regulated genes
        if (length(sig_down) > 5) {
            tryCatch({
                term2gene <- msig[, c("gs_name", "gene_symbol")]
                ora_down <- enricher(gene = sig_down, universe = background,
                                      TERM2GENE = term2gene,
                                      pAdjustMethod = "BH", pvalueCutoff = opt$`padj-threshold`,
                                      minGSSize = opt$`min-size`, maxGSSize = opt$`max-size`)
                if (!is.null(ora_down) && nrow(as.data.frame(ora_down)) > 0) {
                    write.table(as.data.frame(ora_down),
                                file.path(opt$`output-dir`, "ora", paste0(contrast_name, "_", db_name, "_ora_down.tsv")),
                                sep = "\t", quote = FALSE, row.names = FALSE)
                }
            }, error = function(e) cat("    ORA (down) error:", e$message, "\n"))
        }
    }
}

cat("\nPathway enrichment analysis complete.\n")
