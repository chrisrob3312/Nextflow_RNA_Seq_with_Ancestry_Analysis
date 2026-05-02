#!/bin/bash
# ============================================================================
# SETUP DIRECTORIES AND VALIDATE INPUTS
# ============================================================================
# Usage: bash setup_directories.sh
#
# This script (non-SLURM) creates the full output directory tree for the
# cancer RNA-seq pipeline and validates that the samplesheet and key
# reference files exist before any sbatch jobs are submitted.
# ============================================================================

set -euo pipefail

# Source shared configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "============================================================================"
echo "  CANCER RNA-SEQ PIPELINE - DIRECTORY SETUP & INPUT VALIDATION"
echo "============================================================================"
echo ""
echo "Project directory: ${PROJECT_DIR}"
echo "Results directory: ${RESULTS_DIR}"
echo ""

# ============================================================================
# CREATE OUTPUT DIRECTORY TREE
# ============================================================================

echo "[1/3] Creating output directory tree..."

# Top-level directories
mkdir -p "${RESULTS_DIR}"
mkdir -p "${LOGS_DIR}"
mkdir -p "${TEMP_DIR}"

# Pipeline step output directories
mkdir -p "${RESULTS_DIR}/genome_detection"
mkdir -p "${RESULTS_DIR}/qc/samtools"
mkdir -p "${RESULTS_DIR}/qc/rseqc"
mkdir -p "${RESULTS_DIR}/qc/picard"
mkdir -p "${RESULTS_DIR}/qc/multiqc"
mkdir -p "${RESULTS_DIR}/counting/featurecounts"
mkdir -p "${RESULTS_DIR}/counting/matrices"
mkdir -p "${RESULTS_DIR}/counting/normalized"
mkdir -p "${RESULTS_DIR}/ancestry"
mkdir -p "${RESULTS_DIR}/variant_calling"
mkdir -p "${RESULTS_DIR}/hla_typing"
mkdir -p "${RESULTS_DIR}/differential_expression"
mkdir -p "${RESULTS_DIR}/splicing/rmats"
mkdir -p "${RESULTS_DIR}/splicing/leafcutter"
mkdir -p "${RESULTS_DIR}/splicing/spladder"
mkdir -p "${RESULTS_DIR}/splicing/bisbee"
mkdir -p "${RESULTS_DIR}/fusions/fusioncatcher"
mkdir -p "${RESULTS_DIR}/fusions/arriba"
mkdir -p "${RESULTS_DIR}/cnv"
mkdir -p "${RESULTS_DIR}/wgcna"
mkdir -p "${RESULTS_DIR}/pathway_enrichment"
mkdir -p "${RESULTS_DIR}/immune_analysis"
mkdir -p "${RESULTS_DIR}/neoantigen"
mkdir -p "${RESULTS_DIR}/tcr_repertoire"
mkdir -p "${RESULTS_DIR}/sensitivity"
mkdir -p "${RESULTS_DIR}/pharmacogenomics"
mkdir -p "${RESULTS_DIR}/molecular_subtyping"
mkdir -p "${RESULTS_DIR}/validation"
mkdir -p "${RESULTS_DIR}/visualization"

# Log subdirectories
mkdir -p "${LOGS_DIR}/genome_detection"
mkdir -p "${LOGS_DIR}/qc"
mkdir -p "${LOGS_DIR}/counting"
mkdir -p "${LOGS_DIR}/ancestry"
mkdir -p "${LOGS_DIR}/variant_calling"
mkdir -p "${LOGS_DIR}/hla_typing"
mkdir -p "${LOGS_DIR}/differential_expression"
mkdir -p "${LOGS_DIR}/splicing"
mkdir -p "${LOGS_DIR}/fusions"
mkdir -p "${LOGS_DIR}/cnv"
mkdir -p "${LOGS_DIR}/wgcna"
mkdir -p "${LOGS_DIR}/pathway_enrichment"
mkdir -p "${LOGS_DIR}/immune_analysis"
mkdir -p "${LOGS_DIR}/neoantigen"
mkdir -p "${LOGS_DIR}/tcr_repertoire"
mkdir -p "${LOGS_DIR}/sensitivity"
mkdir -p "${LOGS_DIR}/pharmacogenomics"
mkdir -p "${LOGS_DIR}/molecular_subtyping"
mkdir -p "${LOGS_DIR}/validation"
mkdir -p "${LOGS_DIR}/visualization"

echo "  Done. Directory tree created."
echo ""

# ============================================================================
# VALIDATE SAMPLESHEET
# ============================================================================

echo "[2/3] Validating samplesheet..."

ERRORS=0

if [[ ! -f "${SAMPLESHEET}" ]]; then
    echo "  ERROR: Samplesheet not found: ${SAMPLESHEET}"
    ERRORS=$((ERRORS + 1))
else
    echo "  Samplesheet found: ${SAMPLESHEET}"

    # Check header
    HEADER=$(head -n 1 "${SAMPLESHEET}")
    echo "  Header: ${HEADER}"

    # Count samples
    NUM_SAMPLES=$(tail -n +2 "${SAMPLESHEET}" | grep -c . || true)
    echo "  Number of samples: ${NUM_SAMPLES}"

    if [[ "${NUM_SAMPLES}" -eq 0 ]]; then
        echo "  ERROR: Samplesheet contains no sample entries."
        ERRORS=$((ERRORS + 1))
    fi

    # Validate each BAM/BAI file exists
    echo "  Checking BAM/BAI files..."
    MISSING_BAMS=0
    MISSING_BAIS=0

    while IFS=',' read -r sample_id bam_path bai_path rest; do
        [[ -z "${sample_id}" ]] && continue
        if [[ ! -f "${bam_path}" ]]; then
            echo "    WARNING: BAM not found for ${sample_id}: ${bam_path}"
            MISSING_BAMS=$((MISSING_BAMS + 1))
        fi
        if [[ ! -f "${bai_path}" ]]; then
            echo "    WARNING: BAI not found for ${sample_id}: ${bai_path}"
            MISSING_BAIS=$((MISSING_BAIS + 1))
        fi
    done < <(tail -n +2 "${SAMPLESHEET}")

    if [[ "${MISSING_BAMS}" -gt 0 ]]; then
        echo "  WARNING: ${MISSING_BAMS} BAM file(s) not found."
    fi
    if [[ "${MISSING_BAIS}" -gt 0 ]]; then
        echo "  WARNING: ${MISSING_BAIS} BAI file(s) not found."
    fi
fi

echo ""

# ============================================================================
# VALIDATE KEY REFERENCE FILES
# ============================================================================

echo "[3/3] Validating key reference files..."

validate_ref() {
    local path="$1"
    local desc="$2"
    local required="${3:-yes}"

    if [[ -f "${path}" ]]; then
        echo "  OK:      ${desc}"
    elif [[ "${required}" == "yes" ]]; then
        echo "  ERROR:   ${desc} NOT FOUND: ${path}"
        ERRORS=$((ERRORS + 1))
    else
        echo "  SKIPPED: ${desc} (optional, not found: ${path})"
    fi
}

validate_ref "${FASTA}" "Reference genome FASTA"
validate_ref "${FASTA_FAI}" "Reference genome FASTA index"
validate_ref "${GTF}" "Gene annotation GTF"
validate_ref "${GENE_BED}" "Gene BED12 file"
validate_ref "${KNOWN_SNPS}" "Known SNPs VCF (dbSNP)"
validate_ref "${KNOWN_SNPS_TBI}" "Known SNPs VCF index"
validate_ref "${CHAIN_FILE}" "Liftover chain file (hg19->hg38)" "no"
validate_ref "${GRAF_SNP_BED}" "GrafPop SNP file" "no"
validate_ref "${SOMALIER_SITES}" "Somalier sites VCF" "no"
validate_ref "${SOMALIER_ANCESTRY_LABELS}" "Somalier ancestry labels" "no"

echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo "============================================================================"
if [[ "${ERRORS}" -gt 0 ]]; then
    echo "  SETUP COMPLETED WITH ${ERRORS} ERROR(S)."
    echo "  Please fix the errors above before submitting pipeline jobs."
    echo "============================================================================"
    exit 1
else
    echo "  SETUP COMPLETED SUCCESSFULLY."
    echo "  All required files validated. Ready to submit pipeline jobs."
    echo "============================================================================"
    exit 0
fi
