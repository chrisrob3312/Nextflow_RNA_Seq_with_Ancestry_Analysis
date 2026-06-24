#!/bin/bash
#SBATCH --job-name=pharmacogx
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/pharmacogx_%j.out
#SBATCH --error=logs/pharmacogx_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 17 - PharmacoGx Drug Sensitivity Analysis
# ============================================================================
# Uses PharmacoGx PharmacoSets (GDSC2, CTRPv2, gCSI) to train ridge regression
# models for drug sensitivity prediction and compare predictions across ancestry
# groups.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
NORMALIZED_COUNTS="${NORMALIZED_COUNTS:-${RESULTS_DIR}/07_differential_expression/normalized_counts.tsv}"
METADATA="${METADATA:-${SAMPLESHEET}}"
ANCESTRY_PROPORTIONS="${ANCESTRY_PROPORTIONS:-${RESULTS_DIR}/04_ancestry/ancestry_proportions.tsv}"
OUTPUT_DIR="${PHARMACOGX_OUTPUT_DIR:-${RESULTS_DIR}/17_pharmacogenomics/pharmacogx}"
DATASETS="${PHARMACOGX_DATASETS:-GDSC2,CTRPv2,gCSI}"

# --- Validate inputs ---
log_step "PHARMACOGX" "Starting PharmacoGx drug sensitivity analysis"

check_file_exists "${NORMALIZED_COUNTS}" "Normalized count matrix"
check_file_exists "${METADATA}" "Metadata/sample sheet"

if [[ -f "${ANCESTRY_PROPORTIONS}" ]]; then
    log_step "PHARMACOGX" "Ancestry proportions found: ${ANCESTRY_PROPORTIONS}"
else
    log_step "PHARMACOGX" "WARNING: Ancestry proportions not found, ancestry comparisons will be skipped"
    ANCESTRY_PROPORTIONS="none"
fi

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{predictions,comparisons,models,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "PHARMACOGX" "Parameters:"
log_step "PHARMACOGX" "  Normalized counts: ${NORMALIZED_COUNTS}"
log_step "PHARMACOGX" "  Metadata: ${METADATA}"
log_step "PHARMACOGX" "  Ancestry proportions: ${ANCESTRY_PROPORTIONS}"
log_step "PHARMACOGX" "  Output: ${OUTPUT_DIR}"
log_step "PHARMACOGX" "  Datasets: ${DATASETS}"
log_step "PHARMACOGX" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_pharma

# --- Run PharmacoGx analysis R script ---
Rscript "$(dirname "$0")/pharmacogx_analysis.R" \
    --normalized_counts "${NORMALIZED_COUNTS}" \
    --metadata "${METADATA}" \
    --ancestry_proportions "${ANCESTRY_PROPORTIONS}" \
    --output_dir "${OUTPUT_DIR}" \
    --datasets "${DATASETS}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "PHARMACOGX" "PharmacoGx analysis completed successfully"
    log_step "PHARMACOGX" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "PHARMACOGX" "ERROR: PharmacoGx analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "PharmacoGx Analysis Summary"
echo "============================================"
echo "Predictions:"
ls "${OUTPUT_DIR}/predictions/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Ancestry comparisons:"
ls "${OUTPUT_DIR}/comparisons/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Plots:"
ls "${OUTPUT_DIR}/plots/"*.pdf 2>/dev/null | wc -l | xargs -I{} echo "  {} plot files"
echo "============================================"

log_step "PHARMACOGX" "Job complete"
