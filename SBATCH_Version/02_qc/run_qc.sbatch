#!/bin/bash
# ============================================================================
# QC: SAMTOOLS, RSeQC, PICARD, MultiQC
# ============================================================================
# Runs comprehensive QC on STAR-aligned BAM files using SLURM job arrays.
# Each array task processes one sample. After all per-sample tasks complete,
# submit with --dependency=afterok:$JOBID to run MultiQC aggregation.
#
# Submit:
#   # Per-sample QC (array):
#   JOBID=$(sbatch --parsable --array=0-$((N_SAMPLES-1)) run_qc.sbatch)
#
#   # MultiQC aggregation (after all per-sample tasks finish):
#   sbatch --dependency=afterok:${JOBID} run_qc.sbatch --multiqc-only
# ============================================================================

#SBATCH --job-name=rnaseq_qc
#SBATCH --output=../logs/qc/qc_%A_%a.out
#SBATCH --error=../logs/qc/qc_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
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
# module load samtools/1.17
# module load rseqc/5.0.1
# module load picard/3.0.0
# module load multiqc/1.14
# module load python/3.10
# module load R/4.3.0
# conda activate rnaseq_pipeline

# ============================================================================
# DETERMINE MODE: PER-SAMPLE (ARRAY) OR MULTIQC AGGREGATION
# ============================================================================

MULTIQC_ONLY=false
if [[ "${1:-}" == "--multiqc-only" ]]; then
    MULTIQC_ONLY=true
fi

# ============================================================================
# OUTPUT DIRECTORIES
# ============================================================================

SAMTOOLS_DIR="${RESULTS_DIR}/qc/samtools"
RSEQC_DIR="${RESULTS_DIR}/qc/rseqc"
PICARD_DIR="${RESULTS_DIR}/qc/picard"
MULTIQC_DIR="${RESULTS_DIR}/qc/multiqc"

mkdir -p "${SAMTOOLS_DIR}" "${RSEQC_DIR}" "${PICARD_DIR}" "${MULTIQC_DIR}"

# ============================================================================
# MULTIQC AGGREGATION MODE
# ============================================================================

if [[ "${MULTIQC_ONLY}" == "true" ]]; then
    log_step "qc" "Running MultiQC aggregation"

    multiqc \
        --force \
        --outdir "${MULTIQC_DIR}" \
        --title "RNA-Seq QC Report" \
        --comment "Aggregated QC from samtools, RSeQC, and Picard" \
        --module samtools \
        --module rseqc \
        --module picard \
        "${SAMTOOLS_DIR}" \
        "${RSEQC_DIR}" \
        "${PICARD_DIR}"

    log_step "qc" "MultiQC report generated: ${MULTIQC_DIR}/multiqc_report.html"
    echo "MultiQC complete. Report: ${MULTIQC_DIR}/multiqc_report.html"
    exit 0
fi

# ============================================================================
# PER-SAMPLE QC (ARRAY MODE)
# ============================================================================

# Build sample array from samplesheet
mapfile -t SAMPLE_IDS < <(get_sample_ids)

# Get current sample based on array task ID
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID not set. Submit with --array flag or use --multiqc-only."
    echo "Usage: sbatch --array=0-$((${#SAMPLE_IDS[@]}-1)) run_qc.sbatch"
    exit 1
fi

SAMPLE_ID="${SAMPLE_IDS[${SLURM_ARRAY_TASK_ID}]}"
BAM_PATH=$(get_bam_for_sample "${SAMPLE_ID}")

echo "============================================================================"
echo "  QC for sample: ${SAMPLE_ID} (task ${SLURM_ARRAY_TASK_ID})"
echo "  BAM: ${BAM_PATH}"
echo "============================================================================"

# Validate input
if [[ ! -f "${BAM_PATH}" ]]; then
    echo "ERROR: BAM file not found: ${BAM_PATH}"
    exit 1
fi

log_step "qc" "Starting QC for ${SAMPLE_ID}"

# ============================================================================
# SAMTOOLS QC
# ============================================================================

echo ""
echo "--- Samtools flagstat ---"
samtools flagstat \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    "${BAM_PATH}" \
    > "${SAMTOOLS_DIR}/${SAMPLE_ID}.flagstat.txt"

echo "--- Samtools idxstats ---"
samtools idxstats \
    "${BAM_PATH}" \
    > "${SAMTOOLS_DIR}/${SAMPLE_ID}.idxstats.txt"

echo "--- Samtools stats ---"
samtools stats \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    "${BAM_PATH}" \
    > "${SAMTOOLS_DIR}/${SAMPLE_ID}.stats.txt"

echo "  Samtools QC complete."

# ============================================================================
# RSeQC
# ============================================================================

echo ""
echo "--- RSeQC bam_stat.py ---"
bam_stat.py \
    -i "${BAM_PATH}" \
    > "${RSEQC_DIR}/${SAMPLE_ID}.bam_stat.txt" \
    2>&1

echo "--- RSeQC read_distribution.py ---"
read_distribution.py \
    -i "${BAM_PATH}" \
    -r "${GENE_BED}" \
    > "${RSEQC_DIR}/${SAMPLE_ID}.read_distribution.txt" \
    2>&1

echo "--- RSeQC infer_experiment.py ---"
infer_experiment.py \
    -i "${BAM_PATH}" \
    -r "${GENE_BED}" \
    -s 1000000 \
    > "${RSEQC_DIR}/${SAMPLE_ID}.infer_experiment.txt" \
    2>&1

echo "  RSeQC complete."

# ============================================================================
# PICARD CollectRnaSeqMetrics
# ============================================================================

echo ""
echo "--- Picard CollectRnaSeqMetrics ---"

# Determine strandedness for Picard
case "${FC_STRANDEDNESS}" in
    0) PICARD_STRAND="NONE" ;;
    1) PICARD_STRAND="FIRST_READ_TRANSCRIPTION_STRAND" ;;
    2) PICARD_STRAND="SECOND_READ_TRANSCRIPTION_STRAND" ;;
    *) PICARD_STRAND="NONE" ;;
esac

# Generate refFlat from GTF if not present (cached)
REFFLAT="${RESULTS_DIR}/qc/picard/refFlat.txt"
if [[ ! -f "${REFFLAT}" ]]; then
    echo "  Generating refFlat from GTF (first run only)..."
    # Use gtfToGenePred if available, otherwise use a simpler approach
    if command -v gtfToGenePred &> /dev/null; then
        gtfToGenePred -genePredExt -geneNameAsName2 "${GTF}" /dev/stdout | \
            awk 'BEGIN{OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' \
            > "${REFFLAT}"
    else
        echo "  WARNING: gtfToGenePred not found. Using picard CreateSequenceDictionary approach."
        echo "  Ensure refFlat is generated before running this step."
        # Fallback: try using picard's BedToIntervalList approach
        touch "${REFFLAT}"
    fi
fi

# Generate ribosomal intervals if not present (cached)
RRNA_INTERVALS="${RESULTS_DIR}/qc/picard/rRNA.interval_list"
if [[ ! -f "${RRNA_INTERVALS}" ]]; then
    echo "  Generating rRNA interval list (first run only)..."
    # Create sequence dictionary header
    samtools view -H "${BAM_PATH}" | grep "^@SQ" > "${RRNA_INTERVALS}" || true
    samtools view -H "${BAM_PATH}" | grep "^@HD" >> "${RRNA_INTERVALS}" || true
    # Extract rRNA intervals from GTF
    grep -i 'gene_type "rRNA"' "${GTF}" | \
        awk '$3 == "gene" || $3 == "transcript"' | \
        awk 'BEGIN{OFS="\t"} {print $1, $4, $5, $7, "rRNA_"NR}' \
        >> "${RRNA_INTERVALS}" || true
fi

picard CollectRnaSeqMetrics \
    I="${BAM_PATH}" \
    O="${PICARD_DIR}/${SAMPLE_ID}.RNA_Metrics.txt" \
    REF_FLAT="${REFFLAT}" \
    STRAND_SPECIFICITY="${PICARD_STRAND}" \
    RIBOSOMAL_INTERVALS="${RRNA_INTERVALS}" \
    CHART_OUTPUT="${PICARD_DIR}/${SAMPLE_ID}.RNA_Metrics.pdf" \
    VALIDATION_STRINGENCY=LENIENT \
    2>&1 | tee "${PICARD_DIR}/${SAMPLE_ID}.picard.log"

echo "  Picard CollectRnaSeqMetrics complete."

# ============================================================================
# DONE
# ============================================================================

echo ""
echo "============================================================================"
echo "  QC COMPLETE for ${SAMPLE_ID}"
echo "============================================================================"
echo "  Outputs:"
echo "    ${SAMTOOLS_DIR}/${SAMPLE_ID}.flagstat.txt"
echo "    ${SAMTOOLS_DIR}/${SAMPLE_ID}.idxstats.txt"
echo "    ${SAMTOOLS_DIR}/${SAMPLE_ID}.stats.txt"
echo "    ${RSEQC_DIR}/${SAMPLE_ID}.bam_stat.txt"
echo "    ${RSEQC_DIR}/${SAMPLE_ID}.read_distribution.txt"
echo "    ${RSEQC_DIR}/${SAMPLE_ID}.infer_experiment.txt"
echo "    ${PICARD_DIR}/${SAMPLE_ID}.RNA_Metrics.txt"
echo ""

log_step "qc" "QC complete for ${SAMPLE_ID}"
exit 0
