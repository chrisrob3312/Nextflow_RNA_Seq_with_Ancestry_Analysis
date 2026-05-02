#!/bin/bash
# ============================================================================
# FEATURE COUNTING: featureCounts + NORMALIZATION + BATCH CORRECTION
# ============================================================================
# Quantifies gene-level read counts using featureCounts, merges per-sample
# counts into a matrix, and runs normalization/batch correction via R.
#
# Supports two submission modes:
#   1. Array mode (per-sample counting):
#      sbatch --array=0-$((N_SAMPLES-1)) run_featurecounts.sbatch
#
#   2. Merge + normalize mode (after all array tasks complete):
#      sbatch --dependency=afterok:${JOBID} run_featurecounts.sbatch --merge
#
# Or run everything sequentially in a single job (small cohorts):
#      sbatch run_featurecounts.sbatch --all
# ============================================================================

#SBATCH --job-name=featurecounts
#SBATCH --output=../logs/counting/fc_%A_%a.out
#SBATCH --error=../logs/counting/fc_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
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
# MODULE / CONDA SETUP
# ============================================================================

# Uncomment and adjust for your HPC environment:
# module load subread/2.0.6
# module load R/4.3.0
# conda activate rnaseq_pipeline

# ============================================================================
# CONFIGURATION
# ============================================================================

FC_DIR="${RESULTS_DIR}/counting/featurecounts"
MATRIX_DIR="${RESULTS_DIR}/counting/matrices"
NORM_DIR="${RESULTS_DIR}/counting/normalized"

mkdir -p "${FC_DIR}" "${MATRIX_DIR}" "${NORM_DIR}"

# featureCounts parameters from config
STRANDEDNESS="${FC_STRANDEDNESS:-2}"
COUNT_TYPE="${FC_COUNT_TYPE:-exon}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ============================================================================
# DETERMINE MODE
# ============================================================================

MODE="array"
if [[ "${1:-}" == "--merge" ]]; then
    MODE="merge"
elif [[ "${1:-}" == "--all" ]]; then
    MODE="all"
fi

# ============================================================================
# FUNCTION: Run featureCounts for a single sample
# ============================================================================

run_featurecounts_sample() {
    local sample_id="$1"
    local bam_path="$2"

    echo "  Running featureCounts for: ${sample_id}"
    echo "  BAM: ${bam_path}"
    echo "  Strandedness: ${STRANDEDNESS} (0=unstranded, 1=stranded, 2=reversely stranded)"

    if [[ ! -f "${bam_path}" ]]; then
        echo "  ERROR: BAM not found: ${bam_path}"
        return 1
    fi

    featureCounts \
        -a "${GTF}" \
        -o "${FC_DIR}/${sample_id}.featureCounts.txt" \
        -t "${COUNT_TYPE}" \
        -g gene_id \
        -s "${STRANDEDNESS}" \
        -p --countReadPairs \
        -B \
        -C \
        --largestOverlap \
        -T "${THREADS}" \
        --tmpDir "${TEMP_DIR}" \
        "${bam_path}" \
        2>&1 | tee "${FC_DIR}/${sample_id}.featureCounts.log"

    echo "  featureCounts complete for ${sample_id}"
}

# ============================================================================
# FUNCTION: Merge counts into a single matrix
# ============================================================================

merge_counts() {
    echo ""
    echo "============================================================================"
    echo "  MERGING COUNT FILES INTO MATRIX"
    echo "============================================================================"

    local MERGED_COUNTS="${MATRIX_DIR}/raw_counts_matrix.tsv"
    local FIRST=true

    # Get all sample IDs
    mapfile -t SAMPLE_IDS < <(get_sample_ids)

    for sample_id in "${SAMPLE_IDS[@]}"; do
        local fc_file="${FC_DIR}/${sample_id}.featureCounts.txt"

        if [[ ! -f "${fc_file}" ]]; then
            echo "  WARNING: Count file missing for ${sample_id}, skipping."
            continue
        fi

        if [[ "${FIRST}" == "true" ]]; then
            # First sample: take gene_id column + count column
            # featureCounts output: Geneid, Chr, Start, End, Strand, Length, Count
            echo -e "gene_id\t${sample_id}" > "${MERGED_COUNTS}"
            tail -n +3 "${fc_file}" | cut -f1,7 >> "${MERGED_COUNTS}"
            FIRST=false
        else
            # Subsequent samples: paste the count column
            tail -n +3 "${fc_file}" | cut -f7 > "${TEMP_DIR}/${sample_id}.counts_only.tmp"
            paste "${MERGED_COUNTS}" "${TEMP_DIR}/${sample_id}.counts_only.tmp" > "${TEMP_DIR}/merged.tmp"
            mv "${TEMP_DIR}/merged.tmp" "${MERGED_COUNTS}"
            rm -f "${TEMP_DIR}/${sample_id}.counts_only.tmp"

            # Update header
            local CURRENT_HEADER
            CURRENT_HEADER=$(head -n 1 "${MERGED_COUNTS}")
            sed -i "1s/.*/${CURRENT_HEADER}\t${sample_id}/" "${MERGED_COUNTS}"
        fi
    done

    # Also create a gene length file (useful for TPM)
    local GENE_LENGTHS="${MATRIX_DIR}/gene_lengths.tsv"
    local first_sample="${SAMPLE_IDS[0]}"
    local first_fc="${FC_DIR}/${first_sample}.featureCounts.txt"
    if [[ -f "${first_fc}" ]]; then
        echo -e "gene_id\tlength" > "${GENE_LENGTHS}"
        tail -n +3 "${first_fc}" | cut -f1,6 >> "${GENE_LENGTHS}"
    fi

    local N_GENES
    N_GENES=$(tail -n +2 "${MERGED_COUNTS}" | wc -l)
    local N_SAMPLES
    N_SAMPLES=$(head -n 1 "${MERGED_COUNTS}" | awk -F'\t' '{print NF-1}')

    echo "  Merged count matrix: ${MERGED_COUNTS}"
    echo "  Genes: ${N_GENES}"
    echo "  Samples: ${N_SAMPLES}"

    log_step "counting" "Count matrix merged: ${N_GENES} genes x ${N_SAMPLES} samples"
}

# ============================================================================
# FUNCTION: Run normalization and batch correction
# ============================================================================

run_normalization() {
    echo ""
    echo "============================================================================"
    echo "  RUNNING NORMALIZATION AND BATCH CORRECTION"
    echo "============================================================================"

    local RAW_COUNTS="${MATRIX_DIR}/raw_counts_matrix.tsv"
    local R_SCRIPT="${SCRIPT_DIR}/normalize_and_batch_correct.R"

    if [[ ! -f "${RAW_COUNTS}" ]]; then
        echo "ERROR: Raw count matrix not found: ${RAW_COUNTS}"
        exit 1
    fi

    if [[ ! -f "${R_SCRIPT}" ]]; then
        echo "ERROR: R normalization script not found: ${R_SCRIPT}"
        exit 1
    fi

    echo "  Raw counts:          ${RAW_COUNTS}"
    echo "  Metadata:            ${SAMPLESHEET}"
    echo "  Output dir:          ${NORM_DIR}"
    echo "  Normalization:       ${NORMALIZATION_METHOD}"
    echo "  Batch variable:      ${BATCH_VARIABLE}"

    Rscript "${R_SCRIPT}" \
        "${RAW_COUNTS}" \
        "${SAMPLESHEET}" \
        "${NORM_DIR}" \
        "${NORMALIZATION_METHOD}" \
        "${BATCH_VARIABLE}"

    echo "  Normalization and batch correction complete."
    echo "  Outputs in: ${NORM_DIR}/"

    log_step "counting" "Normalization (${NORMALIZATION_METHOD}) and batch correction complete"
}

# ============================================================================
# EXECUTION
# ============================================================================

case "${MODE}" in

    # -------------------------------------------------------------------------
    # ARRAY MODE: one sample per SLURM array task
    # -------------------------------------------------------------------------
    array)
        mapfile -t SAMPLE_IDS < <(get_sample_ids)

        if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
            echo "ERROR: SLURM_ARRAY_TASK_ID not set."
            echo "Submit with: sbatch --array=0-$((${#SAMPLE_IDS[@]}-1)) run_featurecounts.sbatch"
            exit 1
        fi

        SAMPLE_ID="${SAMPLE_IDS[${SLURM_ARRAY_TASK_ID}]}"
        BAM_PATH=$(get_bam_for_sample "${SAMPLE_ID}")

        echo "============================================================================"
        echo "  featureCounts: ${SAMPLE_ID} (array task ${SLURM_ARRAY_TASK_ID})"
        echo "============================================================================"

        log_step "counting" "featureCounts starting for ${SAMPLE_ID}"
        run_featurecounts_sample "${SAMPLE_ID}" "${BAM_PATH}"
        log_step "counting" "featureCounts complete for ${SAMPLE_ID}"
        ;;

    # -------------------------------------------------------------------------
    # MERGE MODE: merge counts + normalize (run after array completes)
    # -------------------------------------------------------------------------
    merge)
        echo "============================================================================"
        echo "  featureCounts: MERGE + NORMALIZE MODE"
        echo "============================================================================"

        merge_counts
        run_normalization
        ;;

    # -------------------------------------------------------------------------
    # ALL MODE: run everything sequentially (small cohorts)
    # -------------------------------------------------------------------------
    all)
        echo "============================================================================"
        echo "  featureCounts: SEQUENTIAL MODE (all samples)"
        echo "============================================================================"

        mapfile -t SAMPLE_IDS < <(get_sample_ids)
        echo "  Processing ${#SAMPLE_IDS[@]} samples sequentially..."
        echo ""

        for sample_id in "${SAMPLE_IDS[@]}"; do
            bam_path=$(get_bam_for_sample "${sample_id}")
            run_featurecounts_sample "${sample_id}" "${bam_path}"
            echo ""
        done

        merge_counts
        run_normalization
        ;;
esac

echo ""
echo "Done."
exit 0
