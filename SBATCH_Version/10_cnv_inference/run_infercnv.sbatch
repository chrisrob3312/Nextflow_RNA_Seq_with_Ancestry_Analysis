#!/bin/bash
# ============================================================================
# COPY NUMBER VARIATION INFERENCE WITH InferCNV
# ============================================================================
# Infers copy number variations from RNA-seq expression data by comparing
# tumor samples against reference (normal) samples using InferCNV with HMM.
# ============================================================================

#SBATCH --job-name=infercnv
#SBATCH --output=logs/infercnv_%j.out
#SBATCH --error=logs/infercnv_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --partition=bigmem
#SBATCH --account=your_account

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load R/4.3.0
# module load jags/4.3.1

log_step "infercnv" "Starting InferCNV analysis"

# --- Create output directories ---
OUTDIR="${RESULTS_DIR}/cnv_inference/infercnv"
mkdir -p "${OUTDIR}"

# --- Create annotations file (tumor vs reference) ---
# Maps each sample to its group: tumor (relapse/no_relapse) or reference (normal)
ANNOTATIONS_FILE="${OUTDIR}/annotations_file.txt"
> "${ANNOTATIONS_FILE}"

REFERENCE_GROUP=""

while IFS=',' read -r sample_id bam bai condition rest; do
    [[ "${sample_id}" == "sample_id" ]] && continue  # skip header

    if [[ "${condition}" == "normal" ]]; then
        echo -e "${sample_id}\treference" >> "${ANNOTATIONS_FILE}"
        REFERENCE_GROUP="reference"
    elif [[ "${condition}" == "relapse" ]]; then
        echo -e "${sample_id}\trelapse" >> "${ANNOTATIONS_FILE}"
    elif [[ "${condition}" == "no_relapse" ]]; then
        echo -e "${sample_id}\tno_relapse" >> "${ANNOTATIONS_FILE}"
    else
        echo -e "${sample_id}\t${condition}" >> "${ANNOTATIONS_FILE}"
    fi
done < "${SAMPLESHEET}"

# If no normal samples, use no_relapse as reference
if [[ -z "${REFERENCE_GROUP}" ]]; then
    log_step "infercnv" "WARNING: No 'normal' samples found. Using 'no_relapse' as reference group."
    REFERENCE_GROUP="no_relapse"
fi

log_step "infercnv" "Annotations file created: ${ANNOTATIONS_FILE}"
log_step "infercnv" "Reference group: ${REFERENCE_GROUP}"

# --- Validate inputs ---
COUNT_MATRIX="${RESULTS_DIR}/counts/raw_counts_matrix.tsv"
GENE_ORDER_FILE="${GENE_BED}"

check_file_exists "${COUNT_MATRIX}" "Raw counts matrix"
check_file_exists "${GENE_ORDER_FILE}" "Gene order/position file"
check_file_exists "${ANNOTATIONS_FILE}" "Cell annotations file"

# --- Run InferCNV via R script ---
log_step "infercnv" "Running InferCNV R script"

Rscript "$(dirname "$0")/run_infercnv.R" \
    --counts_matrix "${COUNT_MATRIX}" \
    --gene_order_file "${GENE_ORDER_FILE}" \
    --annotations_file "${ANNOTATIONS_FILE}" \
    --ref_group "${REFERENCE_GROUP}" \
    --output_dir "${OUTDIR}" \
    --num_threads "${SLURM_CPUS_PER_TASK}"

log_step "infercnv" "InferCNV analysis complete"
log_step "infercnv" "Output directory: ${OUTDIR}"

echo "Done: InferCNV - $(date)"
