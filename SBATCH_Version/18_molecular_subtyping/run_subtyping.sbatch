#!/bin/bash
#SBATCH --job-name=molecular_subtyping
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/subtyping_%j.out
#SBATCH --error=logs/subtyping_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================================
# 18 - Molecular Subtyping (B-ALL Gene Signatures)
# ============================================================================
# Classifies samples into B-ALL molecular subtypes using gene signature
# scoring (GSVA/ssGSEA) across defined subtype signatures.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
NORMALIZED_COUNTS="${NORMALIZED_COUNTS:-${RESULTS_DIR}/07_differential_expression/normalized_counts.tsv}"
METADATA="${METADATA:-${SAMPLESHEET}}"
OUTPUT_DIR="${SUBTYPING_OUTPUT_DIR:-${RESULTS_DIR}/18_molecular_subtyping}"

# --- Validate inputs ---
log_step "SUBTYPING" "Starting molecular subtyping analysis"

check_file_exists "${NORMALIZED_COUNTS}" "Normalized count matrix"
check_file_exists "${METADATA}" "Metadata/sample sheet"

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{scores,assignments,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "SUBTYPING" "Parameters:"
log_step "SUBTYPING" "  Normalized counts: ${NORMALIZED_COUNTS}"
log_step "SUBTYPING" "  Metadata: ${METADATA}"
log_step "SUBTYPING" "  Output: ${OUTPUT_DIR}"
log_step "SUBTYPING" "  CPUs: ${SLURM_CPUS_PER_TASK:-4}"

# --- Load R environment ---
# module load R/4.3.0
# conda activate rnaseq_subtyping

# --- Run molecular subtyping R script ---
Rscript "$(dirname "$0")/molecular_subtyping.R" \
    --normalized_counts "${NORMALIZED_COUNTS}" \
    --metadata "${METADATA}" \
    --output_dir "${OUTPUT_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-4}"

# --- Check completion ---
if [[ $? -eq 0 ]]; then
    log_step "SUBTYPING" "Molecular subtyping completed successfully"
    log_step "SUBTYPING" "Results saved to: ${OUTPUT_DIR}"
else
    log_step "SUBTYPING" "ERROR: Molecular subtyping failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "============================================"
echo "Molecular Subtyping Summary"
echo "============================================"
if [[ -f "${OUTPUT_DIR}/assignments/subtypes.tsv" ]]; then
    echo "Subtype assignments:"
    tail -n +2 "${OUTPUT_DIR}/assignments/subtypes.tsv" | awk -F'\t' '{print $2}' | sort | uniq -c | sort -rn
fi
echo ""
echo "Output files:"
ls "${OUTPUT_DIR}/scores/"*.tsv 2>/dev/null | xargs -I{} echo "  {}"
ls "${OUTPUT_DIR}/assignments/"*.tsv 2>/dev/null | xargs -I{} echo "  {}"
ls "${OUTPUT_DIR}/plots/"*.pdf 2>/dev/null | xargs -I{} echo "  {}"
echo "============================================"

log_step "SUBTYPING" "Job complete"
