#!/bin/bash
#SBATCH --job-name=dream_druggability
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/dream_%j.out
#SBATCH --error=logs/dream_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 17 - DREAM Druggability Analysis
# ============================================================================
# Implements network topology-based druggability scoring: builds co-expression
# networks from DE genes, computes topological metrics (degree, betweenness,
# closeness), scores gene druggability, and identifies drug combination
# candidates from network modules (Louvain community detection).
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
DE_RESULTS_DIR="${DE_RESULTS_DIR:-${RESULTS_DIR}/07_differential_expression}"
NORMALIZED_COUNTS="${NORMALIZED_COUNTS:-${RESULTS_DIR}/07_differential_expression/normalized_counts.tsv}"
OUTPUT_DIR="${DREAM_OUTPUT_DIR:-${RESULTS_DIR}/17_pharmacogenomics/dream}"

# --- Validate inputs ---
log_step "DREAM" "Starting DREAM druggability analysis"

if [[ ! -d "${DE_RESULTS_DIR}" ]]; then
    echo "ERROR: DE results directory not found: ${DE_RESULTS_DIR}"
    exit 1
fi

check_file_exists "${NORMALIZED_COUNTS}" "Normalized count matrix"

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{druggability_scores,combination_candidates,networks,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "DREAM" "Parameters:"
log_step "DREAM" "  DE results: ${DE_RESULTS_DIR}"
log_step "DREAM" "  Normalized counts: ${NORMALIZED_COUNTS}"
log_step "DREAM" "  Output: ${OUTPUT_DIR}"
log_step "DREAM" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_pharma

# --- Run DREAM analysis R script ---
Rscript "$(dirname "$0")/dream_analysis.R" \
    --de_results_dir "${DE_RESULTS_DIR}" \
    --normalized_counts "${NORMALIZED_COUNTS}" \
    --output_dir "${OUTPUT_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "DREAM" "DREAM druggability analysis completed successfully"
    log_step "DREAM" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "DREAM" "ERROR: DREAM druggability analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "DREAM Druggability Summary"
echo "============================================"
echo "Druggability scores:"
ls "${OUTPUT_DIR}/druggability_scores/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Combination candidates:"
ls "${OUTPUT_DIR}/combination_candidates/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Plots:"
ls "${OUTPUT_DIR}/plots/"*.pdf 2>/dev/null | wc -l | xargs -I{} echo "  {} plot files"
echo "============================================"

log_step "DREAM" "Job complete"
