#!/bin/bash
#SBATCH --job-name=wgcna
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/wgcna_%j.out
#SBATCH --error=logs/wgcna_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 11 - Weighted Gene Co-expression Network Analysis (WGCNA)
# ============================================================================
# Builds signed co-expression networks, detects modules, computes eigengenes,
# and tests module-trait correlations including ancestry proportions.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
NORM_COUNTS="${NORM_COUNTS:-${RESULTS_DIR}/07_differential_expression/rdata/normalized_counts.tsv}"
METADATA="${METADATA:-${SAMPLESHEET}}"
ANCESTRY_PROPORTIONS="${ANCESTRY_PROPORTIONS:-${RESULTS_DIR}/04_ancestry/ancestry_proportions.tsv}"
OUTPUT_DIR="${WGCNA_OUTPUT_DIR:-${RESULTS_DIR}/11_wgcna}"
MIN_MODULE_SIZE="${WGCNA_MIN_MODULE_SIZE:-30}"
MERGE_HEIGHT="${WGCNA_MERGE_HEIGHT:-0.25}"
NETWORK_TYPE="${WGCNA_NETWORK_TYPE:-signed}"
TOP_VAR_GENES="${WGCNA_TOP_VAR_GENES:-5000}"

# --- Validate inputs ---
log_step "WGCNA" "Starting WGCNA analysis"

check_file_exists "${NORM_COUNTS}" "Normalized count matrix"
check_file_exists "${METADATA}" "Metadata/sample sheet"

if [[ ! -f "${ANCESTRY_PROPORTIONS}" ]]; then
    echo "WARNING: Ancestry proportions file not found: ${ANCESTRY_PROPORTIONS}"
    echo "         Proceeding without ancestry traits."
    ANCESTRY_PROPORTIONS="none"
fi

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{network,modules,hub_genes,plots,rdata}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "WGCNA" "Parameters:"
log_step "WGCNA" "  Normalized counts: ${NORM_COUNTS}"
log_step "WGCNA" "  Metadata: ${METADATA}"
log_step "WGCNA" "  Ancestry: ${ANCESTRY_PROPORTIONS}"
log_step "WGCNA" "  Output: ${OUTPUT_DIR}"
log_step "WGCNA" "  Network type: ${NETWORK_TYPE}"
log_step "WGCNA" "  Min module size: ${MIN_MODULE_SIZE}"
log_step "WGCNA" "  Merge cut height: ${MERGE_HEIGHT}"
log_step "WGCNA" "  Top variance genes: ${TOP_VAR_GENES}"
log_step "WGCNA" "  CPUs: ${SLURM_CPUS_PER_TASK:-16}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_wgcna

# --- Run WGCNA analysis ---
Rscript "$(dirname "$0")/wgcna_analysis.R" \
    --norm_counts "${NORM_COUNTS}" \
    --metadata "${METADATA}" \
    --ancestry_proportions "${ANCESTRY_PROPORTIONS}" \
    --output_dir "${OUTPUT_DIR}" \
    --min_module_size "${MIN_MODULE_SIZE}" \
    --merge_height "${MERGE_HEIGHT}" \
    --network_type "${NETWORK_TYPE}" \
    --top_var_genes "${TOP_VAR_GENES}" \
    --threads "${SLURM_CPUS_PER_TASK:-16}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "WGCNA" "WGCNA analysis completed successfully"
    log_step "WGCNA" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "WGCNA" "ERROR: WGCNA analysis failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "WGCNA Summary"
echo "============================================"
if [[ -f "${OUTPUT_DIR}/modules/module_assignments.tsv" ]]; then
    n_modules=$(tail -n +2 "${OUTPUT_DIR}/modules/module_assignments.tsv" | cut -f2 | sort -u | wc -l)
    echo "  Modules detected: ${n_modules}"
fi
if [[ -d "${OUTPUT_DIR}/hub_genes" ]]; then
    echo "  Hub gene files:"
    ls "${OUTPUT_DIR}/hub_genes/"*.tsv 2>/dev/null | wc -l | xargs -I{} echo "    {} modules with hub genes"
fi
echo "============================================"

log_step "WGCNA" "Job complete"
