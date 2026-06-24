#!/bin/bash
#SBATCH --job-name=herv_quant
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/herv_%A_%a.out
#SBATCH --error=logs/herv_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-N%10

set -euo pipefail

# ============================================================================
# 14 - HERV/TE Expression Quantification
# ============================================================================
# Array job: per-sample HERV (Human Endogenous Retrovirus) and transposable
# element expression quantification from RNA-seq BAMs.
# Uses Telescope for TE-aware quantification and hervQuant (if reference is
# available) for HERV-specific expression.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
BAM_DIR_INPUT="${BAM_DIR_INPUT:-${BAM_DIR}}"
OUTPUT_DIR="${HERV_OUTPUT_DIR:-${RESULTS_DIR}/neoantigen/herv}"
TELESCOPE_GTF="${TELESCOPE_GTF:-/path/to/refs/retro.hg38.v1.gtf}"
HERVQUANT_REF="${HERVQUANT_REF:-}"

# --- Get sample list ---
mapfile -t SAMPLES < <(get_sample_ids)
NUM_SAMPLES=${#SAMPLES[@]}

SAMPLE_IDX=${SLURM_ARRAY_TASK_ID:-0}
if [[ $SAMPLE_IDX -ge $NUM_SAMPLES ]]; then
    echo "ERROR: Array index ${SAMPLE_IDX} exceeds sample count ${NUM_SAMPLES}"
    exit 1
fi
SAMPLE_ID="${SAMPLES[$SAMPLE_IDX]}"

# --- Validate inputs ---
log_step "HERV" "Processing sample: ${SAMPLE_ID} (index ${SAMPLE_IDX})"

SAMPLE_BAM="${BAM_DIR_INPUT}/${SAMPLE_ID}.bam"
check_file_exists "${SAMPLE_BAM}" "BAM file for ${SAMPLE_ID}"

# --- Create output directories ---
SAMPLE_OUTDIR="${OUTPUT_DIR}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_OUTDIR}"
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "HERV" "Parameters:"
log_step "HERV" "  BAM: ${SAMPLE_BAM}"
log_step "HERV" "  Telescope GTF: ${TELESCOPE_GTF}"
log_step "HERV" "  hervQuant ref: ${HERVQUANT_REF:-not set}"
log_step "HERV" "  Output: ${SAMPLE_OUTDIR}"
log_step "HERV" "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# --- Load software environment ---
# module load telescope/1.0
# module load samtools/1.17
# module load python/3.9
# conda activate telescope

# ============================================================================
# 1. Telescope - TE-aware quantification
# ============================================================================
log_step "HERV" "  Running Telescope for TE quantification..."

if [[ -f "${TELESCOPE_GTF}" ]]; then
    telescope assign \
        "${SAMPLE_BAM}" \
        "${TELESCOPE_GTF}" \
        --attribute locus \
        --no_feature_key __no_feature \
        --ncpu "${SLURM_CPUS_PER_TASK:-8}" \
        --tempdir "${TEMP_DIR}" \
        --outdir "${SAMPLE_OUTDIR}" \
        --exp_tag "${SAMPLE_ID}" \
        --updated_sam

    # Check for output
    TELESCOPE_OUT="${SAMPLE_OUTDIR}/${SAMPLE_ID}-telescope_report.tsv"
    if [[ -f "${TELESCOPE_OUT}" ]]; then
        log_step "HERV" "  Telescope complete: ${TELESCOPE_OUT}"

        # Extract HERV-specific counts
        grep -i "^HERV\|^HML\|^ERV\|^LTR" "${TELESCOPE_OUT}" > \
            "${SAMPLE_OUTDIR}/${SAMPLE_ID}_herv_counts.tsv" 2>/dev/null || true

        HERV_COUNT=$(wc -l < "${SAMPLE_OUTDIR}/${SAMPLE_ID}_herv_counts.tsv" 2>/dev/null || echo 0)
        log_step "HERV" "  HERV elements detected: ${HERV_COUNT}"
    else
        log_step "HERV" "  WARNING: Telescope output not found"
    fi
else
    log_step "HERV" "  WARNING: Telescope GTF not found: ${TELESCOPE_GTF}"
    log_step "HERV" "  Skipping Telescope step"
fi

# ============================================================================
# 2. hervQuant - HERV-specific quantification (if reference available)
# ============================================================================
if [[ -n "${HERVQUANT_REF}" ]] && [[ -d "${HERVQUANT_REF}" ]]; then
    log_step "HERV" "  Running hervQuant..."

    HERVQUANT_OUTDIR="${SAMPLE_OUTDIR}/hervquant"
    mkdir -p "${HERVQUANT_OUTDIR}"

    # hervQuant expects a reference directory with HERV sequences + index
    HERVQUANT_INDEX="${HERVQUANT_REF}/hervquant_index"

    if [[ -f "${HERVQUANT_INDEX}.1.ht2" ]] || [[ -f "${HERVQUANT_INDEX}.1.bt2" ]]; then
        # Extract unmapped / multi-mapped reads for re-alignment to HERV reference
        log_step "HERV" "  Extracting candidate reads..."

        CANDIDATE_R1="${HERVQUANT_OUTDIR}/${SAMPLE_ID}_candidate_R1.fq.gz"
        CANDIDATE_R2="${HERVQUANT_OUTDIR}/${SAMPLE_ID}_candidate_R2.fq.gz"

        samtools view -@ "${SLURM_CPUS_PER_TASK:-8}" -b -f 4 "${SAMPLE_BAM}" | \
            samtools fastq -@ "${SLURM_CPUS_PER_TASK:-8}" \
                -1 "${CANDIDATE_R1}" \
                -2 "${CANDIDATE_R2}" \
                -0 /dev/null \
                -s /dev/null \
                -

        # Align to HERV reference
        log_step "HERV" "  Aligning to HERV reference..."

        HERV_BAM="${HERVQUANT_OUTDIR}/${SAMPLE_ID}_herv_aligned.bam"

        hisat2 \
            -x "${HERVQUANT_INDEX}" \
            -1 "${CANDIDATE_R1}" \
            -2 "${CANDIDATE_R2}" \
            -p "${SLURM_CPUS_PER_TASK:-8}" \
            --no-unal \
            --summary-file "${HERVQUANT_OUTDIR}/${SAMPLE_ID}_herv_alignment_summary.txt" \
            | samtools sort -@ "${SLURM_CPUS_PER_TASK:-8}" -o "${HERV_BAM}"

        samtools index "${HERV_BAM}"

        # Quantify HERV expression (TPM)
        log_step "HERV" "  Quantifying HERV expression..."

        HERV_QUANT_BED="${HERVQUANT_REF}/hervquant_regions.bed"
        if [[ -f "${HERV_QUANT_BED}" ]]; then
            samtools bedcov "${HERV_QUANT_BED}" "${HERV_BAM}" > \
                "${HERVQUANT_OUTDIR}/${SAMPLE_ID}_hervquant_counts.tsv"
        else
            # Fall back to idxstats
            samtools idxstats "${HERV_BAM}" | \
                awk -v s="${SAMPLE_ID}" 'BEGIN{OFS="\t"} $3>0 {print s, $1, $3}' > \
                "${HERVQUANT_OUTDIR}/${SAMPLE_ID}_hervquant_counts.tsv"
        fi

        log_step "HERV" "  hervQuant complete for ${SAMPLE_ID}"

        # Clean up intermediate files
        rm -f "${CANDIDATE_R1}" "${CANDIDATE_R2}"
    else
        log_step "HERV" "  WARNING: hervQuant index not found at ${HERVQUANT_INDEX}"
        log_step "HERV" "  Skipping hervQuant step"
    fi
else
    log_step "HERV" "  hervQuant reference not set or not found, skipping hervQuant"
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "============================================"
echo "HERV/TE Quantification Summary: ${SAMPLE_ID}"
echo "============================================"
echo "Telescope output:"
ls "${SAMPLE_OUTDIR}/"*telescope* 2>/dev/null | head -5
echo "hervQuant output:"
ls "${SAMPLE_OUTDIR}/hervquant/"*hervquant* 2>/dev/null | head -5
echo "============================================"

log_step "HERV" "Sample ${SAMPLE_ID} HERV/TE quantification complete"
