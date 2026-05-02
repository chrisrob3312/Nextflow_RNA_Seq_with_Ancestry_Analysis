#!/bin/bash
#SBATCH --job-name=sensitivity_analysis
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sensitivity_%j.out
#SBATCH --error=logs/sensitivity_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 16 - Sensitivity Analysis
# ============================================================================
# Assesses robustness of differential expression results through:
# 1. Timepoint-stratified DE (diagnostic-only, relapse-only)
# 2. Ancestry threshold sweeps (varying proportion cutoffs)
# 3. Covariate leave-one-out analysis
# 4. Bootstrap DE stability
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
COUNT_MATRIX="${COUNT_MATRIX:-${RESULTS_DIR}/03_counting/raw_counts.tsv}"
METADATA="${METADATA:-${SAMPLESHEET}}"
ANCESTRY_PROPORTIONS="${ANCESTRY_PROPORTIONS:-${RESULTS_DIR}/04_ancestry/ancestry_proportions.tsv}"
DE_RESULTS_DIR="${DE_RESULTS_DIR:-${RESULTS_DIR}/07_differential_expression}"
OUTPUT_DIR="${SENSITIVITY_OUTPUT_DIR:-${RESULTS_DIR}/16_sensitivity}"
BOOT_ITER="${BOOTSTRAP_ITERATIONS:-1000}"
COVARIATES="${DE_COVARIATES:-batch,sex,age,blast_percentage,tumor_purity,timepoint}"

# --- Validate inputs ---
log_step "SENSITIVITY" "Starting sensitivity analysis"

check_file_exists "${COUNT_MATRIX}" "Count matrix"
check_file_exists "${METADATA}" "Metadata/sample sheet"
check_file_exists "${ANCESTRY_PROPORTIONS}" "Ancestry proportions"

if [[ ! -d "${DE_RESULTS_DIR}" ]]; then
    echo "ERROR: DE results directory not found: ${DE_RESULTS_DIR}"
    exit 1
fi

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{timepoint_stratified,ancestry_sweep,leave_one_out,bootstrap,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "SENSITIVITY" "Parameters:"
log_step "SENSITIVITY" "  Count matrix: ${COUNT_MATRIX}"
log_step "SENSITIVITY" "  Metadata: ${METADATA}"
log_step "SENSITIVITY" "  Ancestry proportions: ${ANCESTRY_PROPORTIONS}"
log_step "SENSITIVITY" "  DE results: ${DE_RESULTS_DIR}"
log_step "SENSITIVITY" "  Output: ${OUTPUT_DIR}"
log_step "SENSITIVITY" "  Bootstrap iterations: ${BOOT_ITER}"
log_step "SENSITIVITY" "  Covariates: ${COVARIATES}"
log_step "SENSITIVITY" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_sensitivity

# --- Run sensitivity analysis R script ---
Rscript "$(dirname "$0")/sensitivity_analysis.R" \
    --count_matrix "${COUNT_MATRIX}" \
    --metadata "${METADATA}" \
    --ancestry_proportions "${ANCESTRY_PROPORTIONS}" \
    --de_results_dir "${DE_RESULTS_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --bootstrap_iterations "${BOOT_ITER}" \
    --covariates "${COVARIATES}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "SENSITIVITY" "Sensitivity analysis completed successfully"
    log_step "SENSITIVITY" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "SENSITIVITY" "ERROR: Sensitivity analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "Sensitivity Analysis Summary"
echo "============================================"
echo "Timepoint-stratified results:"
ls "${OUTPUT_DIR}/timepoint_stratified/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Ancestry sweep results:"
ls "${OUTPUT_DIR}/ancestry_sweep/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Leave-one-out results:"
ls "${OUTPUT_DIR}/leave_one_out/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Bootstrap results:"
ls "${OUTPUT_DIR}/bootstrap/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Plots:"
ls "${OUTPUT_DIR}/plots/"*.pdf 2>/dev/null | wc -l | xargs -I{} echo "  {} plot files"
echo "============================================"

log_step "SENSITIVITY" "Job complete"
