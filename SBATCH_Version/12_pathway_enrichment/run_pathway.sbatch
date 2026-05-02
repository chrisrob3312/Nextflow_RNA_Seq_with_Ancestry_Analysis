#!/bin/bash
#SBATCH --job-name=pathway_enrichment
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/pathway_%j.out
#SBATCH --error=logs/pathway_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 12 - Pathway Enrichment Analysis
# ============================================================================
# Runs ORA (clusterProfiler), GSEA (fgsea), and GSVA for pathway-level
# analysis on all DE contrast results. Queries MSigDB collections.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
DE_RESULTS_DIR="${DE_RESULTS_DIR:-${RESULTS_DIR}/07_differential_expression}"
NORM_COUNTS="${NORM_COUNTS:-${RESULTS_DIR}/07_differential_expression/rdata/normalized_counts.tsv}"
METADATA="${METADATA:-${SAMPLESHEET}}"
OUTPUT_DIR="${PATHWAY_OUTPUT_DIR:-${RESULTS_DIR}/12_pathway_enrichment}"
PATHWAY_DATABASES="${PATHWAY_DBS:-GO_BP,GO_MF,KEGG,REACTOME,HALLMARK,IMMUNESIGDB}"
PADJ_THRESHOLD="${DE_PADJ:-0.05}"
LFC_THRESHOLD="${DE_LFC:-0.585}"
SPECIES="${SPECIES:-Homo sapiens}"

# --- Validate inputs ---
log_step "PATHWAY" "Starting pathway enrichment analysis"

check_file_exists "${NORM_COUNTS}" "Normalized count matrix"
check_file_exists "${METADATA}" "Metadata/sample sheet"

if [[ ! -d "${DE_RESULTS_DIR}" ]]; then
    echo "ERROR: DE results directory not found: ${DE_RESULTS_DIR}"
    exit 1
fi

# Check that DE results exist
DE_FILES=$(find "${DE_RESULTS_DIR}/deseq2" -name "*.tsv" 2>/dev/null | wc -l)
if [[ "${DE_FILES}" -eq 0 ]]; then
    echo "ERROR: No DE result files found in ${DE_RESULTS_DIR}/deseq2/"
    exit 1
fi
log_step "PATHWAY" "Found ${DE_FILES} DE result files"

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{ora,gsea,gsva,plots,rdata}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "PATHWAY" "Parameters:"
log_step "PATHWAY" "  DE results dir: ${DE_RESULTS_DIR}"
log_step "PATHWAY" "  Normalized counts: ${NORM_COUNTS}"
log_step "PATHWAY" "  Metadata: ${METADATA}"
log_step "PATHWAY" "  Output: ${OUTPUT_DIR}"
log_step "PATHWAY" "  Pathway databases: ${PATHWAY_DATABASES}"
log_step "PATHWAY" "  padj threshold: ${PADJ_THRESHOLD}"
log_step "PATHWAY" "  LFC threshold: ${LFC_THRESHOLD}"
log_step "PATHWAY" "  Species: ${SPECIES}"
log_step "PATHWAY" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_pathway

# --- Run pathway enrichment R script ---
Rscript "$(dirname "$0")/pathway_enrichment.R" \
    --de_results_dir "${DE_RESULTS_DIR}" \
    --norm_counts "${NORM_COUNTS}" \
    --metadata "${METADATA}" \
    --output_dir "${OUTPUT_DIR}" \
    --pathway_dbs "${PATHWAY_DATABASES}" \
    --padj_threshold "${PADJ_THRESHOLD}" \
    --lfc_threshold "${LFC_THRESHOLD}" \
    --species "${SPECIES}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "PATHWAY" "Pathway enrichment analysis completed successfully"
    log_step "PATHWAY" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "PATHWAY" "ERROR: Pathway enrichment analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "Pathway Enrichment Summary"
echo "============================================"
echo "ORA results:"
ls "${OUTPUT_DIR}/ora/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "GSEA results:"
ls "${OUTPUT_DIR}/gsea/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "GSVA results:"
ls "${OUTPUT_DIR}/gsva/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "============================================"

log_step "PATHWAY" "Job complete"
