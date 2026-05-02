#!/bin/bash
# ============================================================================
# RUN SOMALIER FOR SAMPLE RELATEDNESS AND ANCESTRY VERIFICATION
# ============================================================================
# Extracts genotype information at known informative sites using somalier,
# then runs somalier relate to check sample relatedness and ancestry.
# ============================================================================

#SBATCH --job-name=run_somalier
#SBATCH --output=logs/run_somalier_%j.out
#SBATCH --error=logs/run_somalier_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal
#SBATCH --account=your_account

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load somalier/0.2.16

# --- Setup ---
log_step "somalier" "Starting somalier analysis"

OUTDIR="${RESULTS_DIR}/ancestry/somalier"
EXTRACT_DIR="${OUTDIR}/extracted"
mkdir -p "${OUTDIR}" "${EXTRACT_DIR}"

# --- Validate inputs ---
check_file_exists "${FASTA}" "Reference FASTA"
check_file_exists "${SOMALIER_SITES}" "Somalier sites VCF"

# --- Step 1: somalier extract on each sample ---
log_step "somalier" "Extracting genotypes at informative sites"

mapfile -t SAMPLE_IDS < <(get_sample_ids)

for SAMPLE_ID in "${SAMPLE_IDS[@]}"; do
    BAM_FILE=$(get_bam_for_sample "${SAMPLE_ID}")
    check_file_exists "${BAM_FILE}" "BAM file for ${SAMPLE_ID}"

    log_step "somalier" "Extracting: ${SAMPLE_ID}"

    somalier extract \
        --fasta "${FASTA}" \
        --sites "${SOMALIER_SITES}" \
        --out-dir "${EXTRACT_DIR}" \
        "${BAM_FILE}"
done

log_step "somalier" "Extraction complete for ${#SAMPLE_IDS[@]} samples"

# --- Step 2: somalier relate on all extracted files ---
log_step "somalier" "Running somalier relate"

EXTRACTED_FILES=("${EXTRACT_DIR}"/*.somalier)

if [[ ${#EXTRACTED_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No extracted somalier files found in ${EXTRACT_DIR}" >&2
    exit 1
fi

somalier relate \
    --output-prefix "${OUTDIR}/somalier_relate" \
    "${EXTRACTED_FILES[@]}"

# --- Optionally run ancestry estimation if labels are available ---
if [[ -f "${SOMALIER_ANCESTRY_LABELS}" ]]; then
    log_step "somalier" "Running somalier ancestry estimation"

    somalier relate \
        --output-prefix "${OUTDIR}/somalier_ancestry" \
        --known "${SOMALIER_ANCESTRY_LABELS}" \
        "${EXTRACTED_FILES[@]}"
fi

log_step "somalier" "Somalier analysis complete"
log_step "somalier" "Results in: ${OUTDIR}"
log_step "somalier" "  - somalier_relate.pairs.tsv"
log_step "somalier" "  - somalier_relate.samples.tsv"
log_step "somalier" "  - somalier_relate.html"

echo "Done - $(date)"
