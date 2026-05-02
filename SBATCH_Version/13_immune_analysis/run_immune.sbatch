#!/bin/bash
#SBATCH --job-name=immune_deconv
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/immune_%j.out
#SBATCH --error=logs/immune_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 13 - Immune Deconvolution Analysis
# ============================================================================
# Runs multiple immune deconvolution methods (xCell, MCP-counter, EPIC,
# CIBERSORTx, TIMER) plus ESTIMATE for tumor purity/immune/stromal scoring.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
NORMALIZED_COUNTS="${NORMALIZED_COUNTS:-${RESULTS_DIR}/07_differential_expression/normalized_counts.tsv}"
METADATA="${METADATA:-${SAMPLESHEET}}"
OUTPUT_DIR="${IMMUNE_OUTPUT_DIR:-${RESULTS_DIR}/13_immune_analysis}"
METHODS="${IMMUNE_METHODS:-xcell,mcpcounter,epic,cibersortx,timer}"
CORRECT_TUMOR_PURITY="${CORRECT_TUMOR_PURITY:-TRUE}"

# --- Validate inputs ---
log_step "IMMUNE" "Starting immune deconvolution analysis"

check_file_exists "${NORMALIZED_COUNTS}" "Normalized count matrix"
check_file_exists "${METADATA}" "Metadata/sample sheet"

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{deconvolution,estimate,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "IMMUNE" "Parameters:"
log_step "IMMUNE" "  Normalized counts: ${NORMALIZED_COUNTS}"
log_step "IMMUNE" "  Metadata: ${METADATA}"
log_step "IMMUNE" "  Output: ${OUTPUT_DIR}"
log_step "IMMUNE" "  Methods: ${METHODS}"
log_step "IMMUNE" "  Correct tumor purity: ${CORRECT_TUMOR_PURITY}"
log_step "IMMUNE" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_immune

# --- Run immune deconvolution R script ---
Rscript "$(dirname "$0")/immune_deconvolution.R" \
    --normalized_counts "${NORMALIZED_COUNTS}" \
    --metadata "${METADATA}" \
    --output_dir "${OUTPUT_DIR}" \
    --methods "${METHODS}" \
    --correct_tumor_purity "${CORRECT_TUMOR_PURITY}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "IMMUNE" "Immune deconvolution completed successfully"
    log_step "IMMUNE" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "IMMUNE" "ERROR: Immune deconvolution failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "Immune Deconvolution Summary"
echo "============================================"
echo "Deconvolution results:"
ls "${OUTPUT_DIR}/deconvolution/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "ESTIMATE scores:"
ls "${OUTPUT_DIR}/estimate/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} score files"
echo "Plots:"
ls "${OUTPUT_DIR}/plots/"*.pdf 2>/dev/null | wc -l | xargs -I{} echo "  {} plot files"
echo "============================================"

log_step "IMMUNE" "Job complete"
