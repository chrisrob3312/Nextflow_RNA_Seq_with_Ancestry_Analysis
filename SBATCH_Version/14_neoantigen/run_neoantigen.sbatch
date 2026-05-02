#!/bin/bash
#SBATCH --job-name=neoantigen
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/neoantigen_%A_%a.out
#SBATCH --error=logs/neoantigen_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-N%10

set -euo pipefail

# ============================================================================
# 14 - Neoantigen Prediction (pVACseq + NeoFuse + SNAF)
# ============================================================================
# Array job: each task processes one sample through neoantigen prediction.
# Sources: somatic variants (pVACseq), gene fusions (NeoFuse), and
# aberrant splicing (SNAF).
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
VCF_DIR="${VCF_DIR:-${RESULTS_DIR}/05_variant_calling/annotated_vcfs}"
HLA_DIR="${HLA_DIR:-${RESULTS_DIR}/06_hla_typing}"
FUSION_DIR="${FUSION_DIR:-${RESULTS_DIR}/09_fusions/arriba}"
BAM_DIR_INPUT="${BAM_DIR_INPUT:-${BAM_DIR}}"
OUTPUT_DIR="${NEOANTIGEN_OUTPUT_DIR:-${RESULTS_DIR}/14_neoantigen}"
BINDING_THRESH="${BINDING_THRESHOLD:-500}"

# --- Get sample list ---
mapfile -t SAMPLES < <(get_sample_ids)
NUM_SAMPLES=${#SAMPLES[@]}

# Handle array index
if [[ "${SLURM_ARRAY_TASK_ID:-0}" == "0" ]] && [[ "${1:-}" == "--merge" ]]; then
    # Final merge step (submitted as dependency after array completes)
    echo "=== Running neoantigen merge step ==="
    goto_merge=true
else
    goto_merge=false
    SAMPLE_IDX=${SLURM_ARRAY_TASK_ID:-0}
    if [[ $SAMPLE_IDX -ge $NUM_SAMPLES ]]; then
        echo "ERROR: Array index ${SAMPLE_IDX} exceeds sample count ${NUM_SAMPLES}"
        exit 1
    fi
    SAMPLE_ID="${SAMPLES[$SAMPLE_IDX]}"
fi

# --- Create output directories ---
mkdir -p "${OUTPUT_DIR}"/{pvacseq,neofuse,snaf,merged}
mkdir -p "$(dirname "$0")/logs"

# --- Load software environment ---
# module load pvactools/4.0
# module load neofuse/1.0
# module load samtools/1.17
# conda activate neoantigen

# ============================================================================
# Per-sample neoantigen prediction
# ============================================================================
if [[ "$goto_merge" == "false" ]]; then

    log_step "NEOANTIGEN" "Processing sample: ${SAMPLE_ID} (index ${SAMPLE_IDX})"

    SAMPLE_VCF="${VCF_DIR}/${SAMPLE_ID}.annotated.vcf.gz"
    SAMPLE_BAM="${BAM_DIR_INPUT}/${SAMPLE_ID}.bam"
    SAMPLE_HLA="${HLA_DIR}/${SAMPLE_ID}/hla_alleles.txt"
    SAMPLE_FUSIONS="${FUSION_DIR}/${SAMPLE_ID}/fusions.tsv"
    SAMPLE_OUTDIR="${OUTPUT_DIR}/pvacseq/${SAMPLE_ID}"

    # --- Parse HLA alleles ---
    check_file_exists "${SAMPLE_HLA}" "HLA alleles for ${SAMPLE_ID}"
    HLA_ALLELES=$(cat "${SAMPLE_HLA}" | tr '\n' ',' | sed 's/,$//')
    log_step "NEOANTIGEN" "  HLA alleles: ${HLA_ALLELES}"

    # =========================================================================
    # 1. pVACseq - Somatic variant neoantigens
    # =========================================================================
    if [[ -f "${SAMPLE_VCF}" ]]; then
        log_step "NEOANTIGEN" "  Running pVACseq..."
        mkdir -p "${SAMPLE_OUTDIR}"

        pvacseq run \
            "${SAMPLE_VCF}" \
            "${SAMPLE_ID}" \
            "${HLA_ALLELES}" \
            MHCflurry NetMHCpan NetMHCIIpan \
            "${SAMPLE_OUTDIR}" \
            -e1 8,9,10,11 \
            -e2 15 \
            --binding-threshold "${BINDING_THRESH}" \
            --iedb-install-directory /opt/iedb \
            -t "${SLURM_CPUS_PER_TASK:-8}" \
            --keep-tmp-files

        log_step "NEOANTIGEN" "  pVACseq complete for ${SAMPLE_ID}"
    else
        log_step "NEOANTIGEN" "  WARNING: VCF not found for ${SAMPLE_ID}, skipping pVACseq"
    fi

    # =========================================================================
    # 2. NeoFuse - Fusion-derived neoantigens
    # =========================================================================
    if [[ -f "${SAMPLE_FUSIONS}" ]]; then
        log_step "NEOANTIGEN" "  Running NeoFuse..."
        NEOFUSE_OUTDIR="${OUTPUT_DIR}/neofuse/${SAMPLE_ID}"
        mkdir -p "${NEOFUSE_OUTDIR}"

        NeoFuse \
            -i "${SAMPLE_FUSIONS}" \
            -s "${SAMPLE_ID}" \
            -a "${SAMPLE_HLA}" \
            -o "${NEOFUSE_OUTDIR}" \
            -m 8,9,10,11 \
            -t "${SLURM_CPUS_PER_TASK:-8}" \
            --tool arriba \
            --binding-threshold "${BINDING_THRESH}"

        log_step "NEOANTIGEN" "  NeoFuse complete for ${SAMPLE_ID}"
    else
        log_step "NEOANTIGEN" "  WARNING: Fusion file not found for ${SAMPLE_ID}, skipping NeoFuse"
    fi

    # =========================================================================
    # 3. SNAF - Splicing-derived neoantigens
    # =========================================================================
    if [[ -n "${SNAF_DB:-}" ]] && [[ -d "${SNAF_DB}" ]]; then
        if [[ -f "${SAMPLE_BAM}" ]]; then
            log_step "NEOANTIGEN" "  Running SNAF..."
            SNAF_OUTDIR="${OUTPUT_DIR}/snaf/${SAMPLE_ID}"
            mkdir -p "${SNAF_OUTDIR}"

            # SNAF: identify splicing-derived neoantigens from BAM
            snaf.py \
                --bam "${SAMPLE_BAM}" \
                --gtf "${GTF}" \
                --hla "${SAMPLE_HLA}" \
                --db "${SNAF_DB}" \
                --sample "${SAMPLE_ID}" \
                --output "${SNAF_OUTDIR}" \
                --binding-threshold "${BINDING_THRESH}" \
                --cores "${SLURM_CPUS_PER_TASK:-8}"

            log_step "NEOANTIGEN" "  SNAF complete for ${SAMPLE_ID}"
        else
            log_step "NEOANTIGEN" "  WARNING: BAM not found for ${SAMPLE_ID}, skipping SNAF"
        fi
    else
        log_step "NEOANTIGEN" "  SNAF_DB not set or not found, skipping SNAF"
    fi

    log_step "NEOANTIGEN" "Sample ${SAMPLE_ID} neoantigen prediction complete"
fi

# ============================================================================
# Merge step - Combine all neoantigen sources
# ============================================================================
if [[ "$goto_merge" == "true" ]]; then

    log_step "NEOANTIGEN" "Merging neoantigen predictions across all samples and sources"
    MERGED_DIR="${OUTPUT_DIR}/merged"

    # Create merged summary for each sample
    for SAMPLE_ID in "${SAMPLES[@]}"; do
        MERGED_FILE="${MERGED_DIR}/${SAMPLE_ID}_all_neoantigens.tsv"

        # Header
        echo -e "sample\tsource\tpeptide\thla_allele\taffinity_nm\tgene\tvariant_type" > "${MERGED_FILE}"

        # pVACseq results
        PVAC_RESULT="${OUTPUT_DIR}/pvacseq/${SAMPLE_ID}/MHC_Class_I/${SAMPLE_ID}.filtered.tsv"
        if [[ -f "${PVAC_RESULT}" ]]; then
            tail -n +2 "${PVAC_RESULT}" | awk -F'\t' -v s="${SAMPLE_ID}" \
                '{print s"\tpvacseq\t"$0}' >> "${MERGED_FILE}" || true
        fi

        # NeoFuse results
        NEOFUSE_RESULT="${OUTPUT_DIR}/neofuse/${SAMPLE_ID}/NeoFuse_results.tsv"
        if [[ -f "${NEOFUSE_RESULT}" ]]; then
            tail -n +2 "${NEOFUSE_RESULT}" | awk -F'\t' -v s="${SAMPLE_ID}" \
                '{print s"\tneofuse\t"$0}' >> "${MERGED_FILE}" || true
        fi

        # SNAF results
        SNAF_RESULT="${OUTPUT_DIR}/snaf/${SAMPLE_ID}/snaf_neoantigens.tsv"
        if [[ -f "${SNAF_RESULT}" ]]; then
            tail -n +2 "${SNAF_RESULT}" | awk -F'\t' -v s="${SAMPLE_ID}" \
                '{print s"\tsnaf\t"$0}' >> "${MERGED_FILE}" || true
        fi
    done

    # Create cohort-level summary
    COHORT_SUMMARY="${MERGED_DIR}/neoantigen_summary.tsv"
    echo -e "sample\tpvacseq_count\tneofuse_count\tsnaf_count\ttotal_neoantigens" > "${COHORT_SUMMARY}"

    for SAMPLE_ID in "${SAMPLES[@]}"; do
        PVAC_COUNT=0
        NEOFUSE_COUNT=0
        SNAF_COUNT=0

        PVAC_RESULT="${OUTPUT_DIR}/pvacseq/${SAMPLE_ID}/MHC_Class_I/${SAMPLE_ID}.filtered.tsv"
        [[ -f "${PVAC_RESULT}" ]] && PVAC_COUNT=$(tail -n +2 "${PVAC_RESULT}" | wc -l)

        NEOFUSE_RESULT="${OUTPUT_DIR}/neofuse/${SAMPLE_ID}/NeoFuse_results.tsv"
        [[ -f "${NEOFUSE_RESULT}" ]] && NEOFUSE_COUNT=$(tail -n +2 "${NEOFUSE_RESULT}" | wc -l)

        SNAF_RESULT="${OUTPUT_DIR}/snaf/${SAMPLE_ID}/snaf_neoantigens.tsv"
        [[ -f "${SNAF_RESULT}" ]] && SNAF_COUNT=$(tail -n +2 "${SNAF_RESULT}" | wc -l)

        TOTAL=$((PVAC_COUNT + NEOFUSE_COUNT + SNAF_COUNT))
        echo -e "${SAMPLE_ID}\t${PVAC_COUNT}\t${NEOFUSE_COUNT}\t${SNAF_COUNT}\t${TOTAL}" >> "${COHORT_SUMMARY}"
    done

    log_step "NEOANTIGEN" "Merge complete. Summary: ${COHORT_SUMMARY}"
    echo ""
    echo "============================================"
    echo "Neoantigen Prediction Summary"
    echo "============================================"
    cat "${COHORT_SUMMARY}"
    echo "============================================"
fi

log_step "NEOANTIGEN" "Job complete"
