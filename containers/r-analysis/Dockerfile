FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL description="R environment for cancer RNA-seq analysis pipeline"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'DESeq2', 'edgeR', 'limma', 'sva', \
    'WGCNA', 'dynamicTreeCut', \
    'clusterProfiler', 'fgsea', 'msigdbr', 'enrichplot', 'org.Hs.eg.db', \
    'GSVA', \
    'immunedeconv', 'estimate', \
    'infercnv', \
    'VariantAnnotation', \
    'EnhancedVolcano', \
    'pheatmap', 'RColorBrewer', 'ggrepel', 'patchwork', 'scales', \
    'optparse', 'jsonlite', 'httr' \
    ), update = FALSE, ask = FALSE)"

# CRAN packages
RUN R -e "install.packages(c('oncoPredict', 'pRRophetic'), repos = 'https://cran.r-project.org')" || true

WORKDIR /opt
