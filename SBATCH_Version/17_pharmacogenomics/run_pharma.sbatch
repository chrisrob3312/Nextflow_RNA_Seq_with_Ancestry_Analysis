#!/bin/bash
#SBATCH --job-name=pharmacogenomics
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/pharma_%j.out
#SBATCH --error=logs/pharma_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 17 - Pharmacogenomics Analysis
# ============================================================================
# Predicts drug sensitivity using oncoPredict with GDSC/CCLE training data,
# queries DGIdb for druggable targets, and compares drug sensitivity between
# ancestry/clinical groups.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
NORMALIZED_COUNTS="${NORMALIZED_COUNTS:-${RESULTS_DIR}/07_differential_expression/normalized_counts.tsv}"
DE_RESULTS_DIR="${DE_RESULTS_DIR:-${RESULTS_DIR}/07_differential_expression}"
METADATA="${METADATA:-${SAMPLESHEET}}"
ANCESTRY_PROPORTIONS="${ANCESTRY_PROPORTIONS:-${RESULTS_DIR}/04_ancestry/ancestry_proportions.tsv}"
OUTPUT_DIR="${PHARMA_OUTPUT_DIR:-${RESULTS_DIR}/17_pharmacogenomics}"
DRUG_DB="${DRUG_RESPONSE_DB:-GDSC}"

# --- Validate inputs ---
log_step "PHARMA" "Starting pharmacogenomics analysis"

check_file_exists "${NORMALIZED_COUNTS}" "Normalized count matrix"
check_file_exists "${METADATA}" "Metadata/sample sheet"

if [[ ! -d "${DE_RESULTS_DIR}" ]]; then
    echo "ERROR: DE results directory not found: ${DE_RESULTS_DIR}"
    exit 1
fi

if [[ -f "${ANCESTRY_PROPORTIONS}" ]]; then
    log_step "PHARMA" "Ancestry proportions found: ${ANCESTRY_PROPORTIONS}"
else
    log_step "PHARMA" "WARNING: Ancestry proportions not found, group comparisons will be limited"
    ANCESTRY_PROPORTIONS="none"
fi

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{drug_sensitivity,druggable_targets,group_comparisons,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "PHARMA" "Parameters:"
log_step "PHARMA" "  Normalized counts: ${NORMALIZED_COUNTS}"
log_step "PHARMA" "  DE results: ${DE_RESULTS_DIR}"
log_step "PHARMA" "  Metadata: ${METADATA}"
log_step "PHARMA" "  Ancestry proportions: ${ANCESTRY_PROPORTIONS}"
log_step "PHARMA" "  Output: ${OUTPUT_DIR}"
log_step "PHARMA" "  Drug database: ${DRUG_DB}"
log_step "PHARMA" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_pharma

# --- Run pharmacogenomics R script ---
Rscript "$(dirname "$0")/pharmacogenomics.R" \
    --normalized_counts "${NORMALIZED_COUNTS}" \
    --de_results_dir "${DE_RESULTS_DIR}" \
    --metadata "${METADATA}" \
    --ancestry_proportions "${ANCESTRY_PROPORTIONS}" \
    --output_dir "${OUTPUT_DIR}" \
    --drug_db "${DRUG_DB}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "PHARMA" "Pharmacogenomics analysis completed successfully"
    log_step "PHARMA" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "PHARMA" "ERROR: Pharmacogenomics analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "Pharmacogenomics Summary"
echo "============================================"
echo "Drug sensitivity predictions:"
ls "${OUTPUT_DIR}/drug_sensitivity/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Druggable targets:"
ls "${OUTPUT_DIR}/druggable_targets/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Group comparisons:"
ls "${OUTPUT_DIR}/group_comparisons/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Plots:"
ls "${OUTPUT_DIR}/plots/"*.pdf 2>/dev/null | wc -l | xargs -I{} echo "  {} plot files"
echo "============================================"

log_step "PHARMA" "Job complete"
