#!/bin/bash
# ============================================================================
# PIPELINE VALIDATION CHECKS
# ============================================================================
# Runs cross-method concordance, sex checks, outlier detection, and genomic
# inflation analysis to validate pipeline outputs before final reporting.
#
# Submit: sbatch run_validation.sbatch
# ============================================================================

#SBATCH --job-name=validation
#SBATCH --output=../logs/validation/validation_%j.out
#SBATCH --error=../logs/validation/validation_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL,END
# #SBATCH --account=${SLURM_ACCOUNT}

set -euo pipefail

# ============================================================================
# SOURCE CONFIG
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# ============================================================================
# PARAMETERS
# ============================================================================

DESEQ2_DIR="${DESEQ2_DIR:-${RESULTS_DIR}/07_differential_expression/deseq2}"
LIMMA_DIR="${LIMMA_DIR:-${RESULTS_DIR}/07_differential_expression/limma_voom}"
COUNT_MATRIX="${COUNT_MATRIX:-${RESULTS_DIR}/03_counting/raw_counts.tsv}"
NORMALIZED_COUNTS="${NORMALIZED_COUNTS:-${RESULTS_DIR}/03_counting/normalized/vst_counts.tsv}"
METADATA="${METADATA:-${SAMPLESHEET}}"
OUTPUT_DIR="${VALIDATION_OUTPUT_DIR:-${RESULTS_DIR}/19_validation}"
OUTLIER_SD="${OUTLIER_SD_THRESHOLD:-3}"
MIN_CORR="${MIN_SAMPLE_CORRELATION:-0.8}"

# ============================================================================
# VALIDATE INPUTS
# ============================================================================

log_step "VALIDATION" "Starting pipeline validation checks"

check_file_exists "${COUNT_MATRIX}" "Raw count matrix"
check_file_exists "${NORMALIZED_COUNTS}" "Normalized count matrix"
check_file_exists "${METADATA}" "Metadata/samplesheet"

if [[ ! -d "${DESEQ2_DIR}" ]]; then
    echo "WARNING: DESeq2 results directory not found: ${DESEQ2_DIR}"
    echo "         DE concordance check will be skipped."
fi

if [[ ! -d "${LIMMA_DIR}" ]]; then
    echo "WARNING: limma-voom results directory not found: ${LIMMA_DIR}"
    echo "         DE concordance check will be skipped."
fi

# ============================================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================================

mkdir -p "${OUTPUT_DIR}"/{plots,reports}
mkdir -p "$(dirname "$0")/logs"

# ============================================================================
# LOG PARAMETERS
# ============================================================================

log_step "VALIDATION" "Parameters:"
log_step "VALIDATION" "  DESeq2 dir: ${DESEQ2_DIR}"
log_step "VALIDATION" "  limma dir: ${LIMMA_DIR}"
log_step "VALIDATION" "  Count matrix: ${COUNT_MATRIX}"
log_step "VALIDATION" "  Normalized counts: ${NORMALIZED_COUNTS}"
log_step "VALIDATION" "  Metadata: ${METADATA}"
log_step "VALIDATION" "  Output dir: ${OUTPUT_DIR}"
log_step "VALIDATION" "  Outlier SD threshold: ${OUTLIER_SD}"
log_step "VALIDATION" "  Min sample correlation: ${MIN_CORR}"
log_step "VALIDATION" "  CPUs: ${SLURM_CPUS_PER_TASK:-4}"

# ============================================================================
# LOAD R ENVIRONMENT
# ============================================================================

# Uncomment whichever module/conda system your cluster uses:
# module load R/4.3.0
# conda activate rnaseq_validation

# ============================================================================
# RUN VALIDATION
# ============================================================================

Rscript "${SCRIPT_DIR}/validation.R" \
    --deseq2_results_dir "${DESEQ2_DIR}" \
    --limma_results_dir "${LIMMA_DIR}" \
    --count_matrix "${COUNT_MATRIX}" \
    --normalized_counts "${NORMALIZED_COUNTS}" \
    --metadata "${METADATA}" \
    --output_dir "${OUTPUT_DIR}" \
    --outlier_sd_threshold "${OUTLIER_SD}" \
    --min_sample_correlation "${MIN_CORR}" \
    --threads "${SLURM_CPUS_PER_TASK:-4}"

# ============================================================================
# CHECK COMPLETION
# ============================================================================

if [[ $? -eq 0 ]]; then
    log_step "VALIDATION" "Validation checks completed successfully"
    log_step "VALIDATION" "Report: ${OUTPUT_DIR}/reports/validation_report.md"
else
    log_step "VALIDATION" "ERROR: Validation checks failed"
    exit 1
fi

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "============================================"
echo "Validation Summary"
echo "============================================"
if [[ -f "${OUTPUT_DIR}/reports/validation_report.md" ]]; then
    echo "Report generated: ${OUTPUT_DIR}/reports/validation_report.md"
    # Print any critical warnings
    grep -c "FAIL\|WARNING\|MISMATCH" "${OUTPUT_DIR}/reports/validation_report.md" 2>/dev/null | \
        xargs -I{} echo "  {} potential issues flagged"
fi
echo "============================================"

log_step "VALIDATION" "Job complete"
