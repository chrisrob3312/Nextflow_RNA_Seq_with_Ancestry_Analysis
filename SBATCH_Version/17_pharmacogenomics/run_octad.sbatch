#!/bin/bash
#SBATCH --job-name=octad
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/octad_%j.out
#SBATCH --error=logs/octad_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 17 - OCTAD Drug Reversal Analysis
# ============================================================================
# Uses the OCTAD Bioconductor package to identify drugs that reverse disease
# gene expression signatures via the LINCS L1000 connectivity map (sRGES).
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
DE_RESULTS_DIR="${DE_RESULTS_DIR:-${RESULTS_DIR}/07_differential_expression}"
OUTPUT_DIR="${OCTAD_OUTPUT_DIR:-${RESULTS_DIR}/17_pharmacogenomics/octad}"

# --- Validate inputs ---
log_step "OCTAD" "Starting OCTAD drug reversal analysis"

if [[ ! -d "${DE_RESULTS_DIR}" ]]; then
    echo "ERROR: DE results directory not found: ${DE_RESULTS_DIR}"
    exit 1
fi

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{drug_candidates,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "OCTAD" "Parameters:"
log_step "OCTAD" "  DE results: ${DE_RESULTS_DIR}"
log_step "OCTAD" "  Output: ${OUTPUT_DIR}"
log_step "OCTAD" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_pharma

# --- Run OCTAD analysis R script ---
Rscript "$(dirname "$0")/octad_analysis.R" \
    --de_results_dir "${DE_RESULTS_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "OCTAD" "OCTAD analysis completed successfully"
    log_step "OCTAD" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "OCTAD" "ERROR: OCTAD analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "OCTAD Analysis Summary"
echo "============================================"
echo "Drug candidate files:"
ls "${OUTPUT_DIR}/drug_candidates/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Plots:"
ls "${OUTPUT_DIR}/plots/"*.pdf 2>/dev/null | wc -l | xargs -I{} echo "  {} plot files"
echo "============================================"

log_step "OCTAD" "Job complete"
