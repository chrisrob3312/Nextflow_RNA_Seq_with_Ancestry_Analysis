#!/bin/bash
#SBATCH --job-name=trust4_tcr
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/trust4_%A_%a.out
#SBATCH --error=logs/trust4_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-N%10

set -euo pipefail

# ============================================================================
# 15 - TCR/BCR Repertoire Reconstruction (TRUST4)
# ============================================================================
# Array job: each task reconstructs immune receptor repertoire from RNA-seq BAM.
# Final merge step computes diversity metrics across the cohort.
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
BAM_DIR_INPUT="${BAM_DIR_INPUT:-${BAM_DIR}}"
OUTPUT_DIR="${TCR_OUTPUT_DIR:-${RESULTS_DIR}/15_tcr_repertoire}"
TRUST4_REF="${TRUST4_REF:-/path/to/refs/human_IMGT+C.fa}"

# --- Get sample list ---
mapfile -t SAMPLES < <(get_sample_ids)
NUM_SAMPLES=${#SAMPLES[@]}

# Handle array index vs merge step
if [[ "${SLURM_ARRAY_TASK_ID:-0}" == "0" ]] && [[ "${1:-}" == "--merge" ]]; then
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
mkdir -p "${OUTPUT_DIR}"/{per_sample,diversity,plots}
mkdir -p "$(dirname "$0")/logs"

# --- Load software environment ---
# module load trust4/1.0.7
# module load samtools/1.17
# conda activate tcr_repertoire

# ============================================================================
# Per-sample TRUST4 repertoire reconstruction
# ============================================================================
if [[ "$goto_merge" == "false" ]]; then

    log_step "TRUST4" "Processing sample: ${SAMPLE_ID} (index ${SAMPLE_IDX})"

    SAMPLE_BAM="${BAM_DIR_INPUT}/${SAMPLE_ID}.bam"
    SAMPLE_OUTDIR="${OUTPUT_DIR}/per_sample/${SAMPLE_ID}"
    mkdir -p "${SAMPLE_OUTDIR}"

    # --- Validate inputs ---
    check_file_exists "${SAMPLE_BAM}" "BAM file for ${SAMPLE_ID}"
    check_file_exists "${TRUST4_REF}" "TRUST4 IMGT reference"
    check_file_exists "${FASTA}" "Reference genome FASTA"

    # --- Run TRUST4 ---
    log_step "TRUST4" "  Running TRUST4 on ${SAMPLE_ID}..."

    run-trust4 \
        -b "${SAMPLE_BAM}" \
        -f "${FASTA}" \
        --ref "${TRUST4_REF}" \
        -o "${SAMPLE_OUTDIR}/${SAMPLE_ID}" \
        -t "${SLURM_CPUS_PER_TASK:-8}"

    # --- Verify outputs ---
    TRUST4_REPORT="${SAMPLE_OUTDIR}/${SAMPLE_ID}_report.tsv"
    if [[ -f "${TRUST4_REPORT}" ]]; then
        TOTAL_CLONES=$(tail -n +2 "${TRUST4_REPORT}" | wc -l)
        log_step "TRUST4" "  TRUST4 complete: ${TOTAL_CLONES} clonotypes identified"
    else
        log_step "TRUST4" "  WARNING: TRUST4 report not generated for ${SAMPLE_ID}"
    fi

    log_step "TRUST4" "Sample ${SAMPLE_ID} complete"
fi

# ============================================================================
# Merge step - Compute diversity metrics across cohort
# ============================================================================
if [[ "$goto_merge" == "true" ]]; then

    log_step "TRUST4" "Computing diversity metrics across all samples"
    DIVERSITY_FILE="${OUTPUT_DIR}/diversity/tcr_diversity_metrics.tsv"

    # Header for diversity metrics
    echo -e "sample_id\tchain\ttotal_clonotypes\tunique_CDR3\tshannon_diversity\tsimpson_diversity\tclonality\tmax_clone_fraction" > "${DIVERSITY_FILE}"

    for SAMPLE_ID in "${SAMPLES[@]}"; do
        TRUST4_REPORT="${OUTPUT_DIR}/per_sample/${SAMPLE_ID}/${SAMPLE_ID}_report.tsv"

        if [[ ! -f "${TRUST4_REPORT}" ]]; then
            log_step "TRUST4" "  WARNING: No report for ${SAMPLE_ID}, skipping"
            continue
        fi

        # Compute diversity metrics per chain (TRA, TRB, IGH, IGK, IGL)
        for CHAIN in TRA TRB IGH IGK IGL; do
            # Extract chain-specific clonotypes
            CHAIN_DATA=$(awk -F'\t' -v chain="${CHAIN}" '$1 ~ chain {print $0}' "${TRUST4_REPORT}")

            if [[ -z "${CHAIN_DATA}" ]]; then
                continue
            fi

            TOTAL_CLONES=$(echo "${CHAIN_DATA}" | wc -l)
            UNIQUE_CDR3=$(echo "${CHAIN_DATA}" | awk -F'\t' '{print $3}' | sort -u | wc -l)

            # Compute Shannon and Simpson diversity using awk
            # Extract abundance/frequency column and compute metrics
            DIVERSITY_STATS=$(echo "${CHAIN_DATA}" | awk -F'\t' '
            BEGIN {
                total = 0; shannon = 0; simpson = 0; max_frac = 0;
            }
            {
                counts[NR] = $2 + 0;  # abundance column
                total += counts[NR];
            }
            END {
                if (total == 0) {
                    print "0\t0\t0\t0";
                    exit;
                }
                for (i = 1; i <= NR; i++) {
                    p = counts[i] / total;
                    if (p > 0) {
                        shannon -= p * log(p) / log(2);
                        simpson += p * p;
                    }
                    if (p > max_frac) max_frac = p;
                }
                max_shannon = (NR > 1) ? log(NR) / log(2) : 1;
                clonality = (max_shannon > 0) ? 1 - (shannon / max_shannon) : 0;
                printf "%.4f\t%.4f\t%.4f\t%.4f\n", shannon, 1 - simpson, clonality, max_frac;
            }')

            echo -e "${SAMPLE_ID}\t${CHAIN}\t${TOTAL_CLONES}\t${UNIQUE_CDR3}\t${DIVERSITY_STATS}" >> "${DIVERSITY_FILE}"
        done
    done

    # Create summary across samples
    SUMMARY_FILE="${OUTPUT_DIR}/diversity/repertoire_summary.tsv"
    echo -e "sample_id\ttotal_TCR_clones\ttotal_BCR_clones\tTRB_diversity\tIGH_diversity" > "${SUMMARY_FILE}"

    for SAMPLE_ID in "${SAMPLES[@]}"; do
        TCR_CLONES=$(awk -F'\t' -v s="${SAMPLE_ID}" '$1==s && ($2=="TRA" || $2=="TRB") {sum+=$3} END{print sum+0}' "${DIVERSITY_FILE}")
        BCR_CLONES=$(awk -F'\t' -v s="${SAMPLE_ID}" '$1==s && ($2=="IGH" || $2=="IGK" || $2=="IGL") {sum+=$3} END{print sum+0}' "${DIVERSITY_FILE}")
        TRB_DIV=$(awk -F'\t' -v s="${SAMPLE_ID}" '$1==s && $2=="TRB" {print $5}' "${DIVERSITY_FILE}" | head -1)
        IGH_DIV=$(awk -F'\t' -v s="${SAMPLE_ID}" '$1==s && $2=="IGH" {print $5}' "${DIVERSITY_FILE}" | head -1)
        echo -e "${SAMPLE_ID}\t${TCR_CLONES}\t${BCR_CLONES}\t${TRB_DIV:-NA}\t${IGH_DIV:-NA}" >> "${SUMMARY_FILE}"
    done

    log_step "TRUST4" "Diversity metrics saved to: ${DIVERSITY_FILE}"
    log_step "TRUST4" "Repertoire summary saved to: ${SUMMARY_FILE}"

    echo ""
    echo "============================================"
    echo "TCR/BCR Repertoire Summary"
    echo "============================================"
    echo "Samples processed: ${NUM_SAMPLES}"
    echo "Diversity metrics: ${DIVERSITY_FILE}"
    echo "Summary: ${SUMMARY_FILE}"
    echo "============================================"
fi

log_step "TRUST4" "Job complete"
