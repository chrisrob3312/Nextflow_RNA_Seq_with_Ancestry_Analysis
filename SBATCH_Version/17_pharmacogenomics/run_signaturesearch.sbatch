#!/bin/bash
#SBATCH --job-name=signaturesearch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/signaturesearch_%j.out
#SBATCH --error=logs/signaturesearch_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 17 - signatureSearch Drug Connectivity Analysis
# ============================================================================
# Uses signatureSearch + signatureSearchData to query LINCS L1000 with disease
# gene signatures (CMAP method), then performs Drug Set Enrichment Analysis
# (DSEA) by mechanism of action.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
DE_RESULTS_DIR="${DE_RESULTS_DIR:-${RESULTS_DIR}/07_differential_expression}"
OUTPUT_DIR="${SIGSEARCH_OUTPUT_DIR:-${RESULTS_DIR}/17_pharmacogenomics/signaturesearch}"

# --- Validate inputs ---
log_step "SIGSEARCH" "Starting signatureSearch analysis"

if [[ ! -d "${DE_RESULTS_DIR}" ]]; then
    echo "ERROR: DE results directory not found: ${DE_RESULTS_DIR}"
    exit 1
fi

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{drug_rankings,dsea,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "SIGSEARCH" "Parameters:"
log_step "SIGSEARCH" "  DE results: ${DE_RESULTS_DIR}"
log_step "SIGSEARCH" "  Output: ${OUTPUT_DIR}"
log_step "SIGSEARCH" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_pharma

# --- Run signatureSearch analysis R script ---
Rscript "$(dirname "$0")/signaturesearch_analysis.R" \
    --de_results_dir "${DE_RESULTS_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "SIGSEARCH" "signatureSearch analysis completed successfully"
    log_step "SIGSEARCH" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "SIGSEARCH" "ERROR: signatureSearch analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "signatureSearch Analysis Summary"
echo "============================================"
echo "Drug rankings:"
ls "${OUTPUT_DIR}/drug_rankings/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "DSEA results:"
ls "${OUTPUT_DIR}/dsea/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Plots:"
ls "${OUTPUT_DIR}/plots/"*.pdf 2>/dev/null | wc -l | xargs -I{} echo "  {} plot files"
echo "============================================"

log_step "SIGSEARCH" "Job complete"
