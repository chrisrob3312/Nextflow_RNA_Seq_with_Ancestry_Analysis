#!/bin/bash
#SBATCH --job-name=neoantigen_burden
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/neoantigen_burden_%j.out
#SBATCH --error=logs/neoantigen_burden_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 14 - Neoantigen Burden Analysis
# ============================================================================
# Calculates per-sample neoantigen burden by source (SNV, fusion, splicing,
# HERV) and compares across clinical and ancestry groups. Generates
# statistical comparisons and visualizations.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
MERGED_NEOANTIGENS="${MERGED_NEOANTIGENS:-${RESULTS_DIR}/14_neoantigen/merged/neoantigen_summary.tsv}"
METADATA="${METADATA:-${SAMPLESHEET}}"
ANCESTRY_PROPORTIONS="${ANCESTRY_PROPORTIONS:-${RESULTS_DIR}/04_ancestry/ancestry_proportions.tsv}"
OUTPUT_DIR="${BURDEN_OUTPUT_DIR:-${RESULTS_DIR}/14_neoantigen/burden_analysis}"

# --- Validate inputs ---
log_step "NEO_BURDEN" "Starting neoantigen burden analysis"

check_file_exists "${MERGED_NEOANTIGENS}" "Merged neoantigen summary"
check_file_exists "${METADATA}" "Metadata/sample sheet"

if [[ -f "${ANCESTRY_PROPORTIONS}" ]]; then
    log_step "NEO_BURDEN" "Ancestry proportions found: ${ANCESTRY_PROPORTIONS}"
else
    log_step "NEO_BURDEN" "WARNING: Ancestry proportions not found, ancestry comparisons will be skipped"
    ANCESTRY_PROPORTIONS="none"
fi

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{comparisons,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "NEO_BURDEN" "Parameters:"
log_step "NEO_BURDEN" "  Merged neoantigens: ${MERGED_NEOANTIGENS}"
log_step "NEO_BURDEN" "  Metadata: ${METADATA}"
log_step "NEO_BURDEN" "  Ancestry proportions: ${ANCESTRY_PROPORTIONS}"
log_step "NEO_BURDEN" "  Output: ${OUTPUT_DIR}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_neoantigen

# --- Run neoantigen burden R script ---
Rscript "$(dirname "$0")/neoantigen_burden.R" \
    --merged_neoantigens "${MERGED_NEOANTIGENS}" \
    --metadata "${METADATA}" \
    --ancestry_proportions "${ANCESTRY_PROPORTIONS}" \
    --output_dir "${OUTPUT_DIR}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "NEO_BURDEN" "Neoantigen burden analysis completed successfully"
    log_step "NEO_BURDEN" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "NEO_BURDEN" "ERROR: Neoantigen burden analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "Neoantigen Burden Analysis Summary"
echo "============================================"
echo "Comparisons:"
ls "${OUTPUT_DIR}/comparisons/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "  {} result files"
echo "Plots:"
ls "${OUTPUT_DIR}/plots/"*.pdf 2>/dev/null | wc -l | xargs -I{} echo "  {} plot files"
echo "Burden table:"
ls "${OUTPUT_DIR}/"*burden*.tsv 2>/dev/null | head -3
echo "============================================"

log_step "NEO_BURDEN" "Job complete"
