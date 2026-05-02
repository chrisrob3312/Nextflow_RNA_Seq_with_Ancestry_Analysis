#!/bin/bash
#SBATCH --job-name=differential_expression
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/de_%j.out
#SBATCH --error=logs/de_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 07 - Differential Expression Analysis (DESeq2 + limma-voom)
# ============================================================================
# Runs both DESeq2 and limma-voom on all contrasts defined in a JSON file.
# Incorporates ancestry proportions as covariates.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters (override via environment or defaults from config) ---
COUNT_MATRIX="${COUNT_MATRIX:-${RESULTS_DIR}/03_counting/raw_counts.tsv}"
METADATA="${METADATA:-${SAMPLESHEET}}"
ANCESTRY_PROPORTIONS="${ANCESTRY_PROPORTIONS:-${RESULTS_DIR}/04_ancestry/ancestry_proportions.tsv}"
CONTRASTS_JSON="${CONTRASTS_JSON:-${PROJECT_DIR}/contrasts.json}"
COVARIATES="${DE_COVARIATES:-batch,sex,age}"
OUTPUT_DIR="${DE_OUTPUT_DIR:-${RESULTS_DIR}/07_differential_expression}"
PADJ_THRESHOLD="${DE_PADJ:-0.05}"
LFC_THRESHOLD="${DE_LFC:-0.585}"

# --- Validate inputs ---
log_step "DE" "Starting differential expression analysis"

check_file_exists "${COUNT_MATRIX}" "Count matrix"
check_file_exists "${METADATA}" "Metadata/sample sheet"
check_file_exists "${CONTRASTS_JSON}" "Contrasts JSON"

if [[ ! -f "${ANCESTRY_PROPORTIONS}" ]]; then
    echo "WARNING: Ancestry proportions file not found: ${ANCESTRY_PROPORTIONS}"
    echo "         Proceeding without ancestry covariates."
    ANCESTRY_PROPORTIONS="none"
fi

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{deseq2,limma_voom,plots,rdata}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "DE" "Parameters:"
log_step "DE" "  Count matrix: ${COUNT_MATRIX}"
log_step "DE" "  Metadata: ${METADATA}"
log_step "DE" "  Ancestry: ${ANCESTRY_PROPORTIONS}"
log_step "DE" "  Contrasts: ${CONTRASTS_JSON}"
log_step "DE" "  Covariates: ${COVARIATES}"
log_step "DE" "  Output: ${OUTPUT_DIR}"
log_step "DE" "  padj threshold: ${PADJ_THRESHOLD}"
log_step "DE" "  LFC threshold: ${LFC_THRESHOLD}"
log_step "DE" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load R environment ---
# Uncomment whichever module/conda system your cluster uses:
# module load R/4.3.0
# conda activate rnaseq_de

# --- Run the DE analysis R script ---
Rscript "$(dirname "$0")/run_de_analysis.R" \
    --count_matrix "${COUNT_MATRIX}" \
    --metadata "${METADATA}" \
    --ancestry_proportions "${ANCESTRY_PROPORTIONS}" \
    --contrasts_json "${CONTRASTS_JSON}" \
    --covariates "${COVARIATES}" \
    --output_dir "${OUTPUT_DIR}" \
    --padj_threshold "${PADJ_THRESHOLD}" \
    --lfc_threshold "${LFC_THRESHOLD}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "DE" "Differential expression analysis completed successfully"
    log_step "DE" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "DE" "ERROR: Differential expression analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "DE Analysis Summary"
echo "============================================"
if [[ -d "${OUTPUT_DIR}/deseq2" ]]; then
    echo "DESeq2 results:"
    ls "${OUTPUT_DIR}/deseq2/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} contrast result files"
fi
if [[ -d "${OUTPUT_DIR}/limma_voom" ]]; then
    echo "limma-voom results:"
    ls "${OUTPUT_DIR}/limma_voom/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} contrast result files"
fi
echo "============================================"

log_step "DE" "Job complete"
