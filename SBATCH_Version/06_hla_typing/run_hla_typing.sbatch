#!/bin/bash
# ============================================================================
# HLA TYPING WITH arcasHLA
# ============================================================================
# SLURM array job: per-sample HLA typing using arcasHLA.
# Extracts HLA reads from STAR-aligned BAMs and genotypes class I + II alleles.
# After all array tasks complete, run this script again with --merge flag
# to consolidate all genotypes.
# ============================================================================

#SBATCH --job-name=hla_typing
#SBATCH --output=logs/hla_typing_%A_%a.out
#SBATCH --error=logs/hla_typing_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --account=your_account
#SBATCH --array=0-N

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load arcasHLA/0.5.0
# module load samtools/1.18

# --- Create output directories ---
OUTDIR="${RESULTS_DIR}/hla_typing"
EXTRACT_DIR="${OUTDIR}/extracted"
GENOTYPE_DIR="${OUTDIR}/genotypes"
mkdir -p "${OUTDIR}" "${EXTRACT_DIR}" "${GENOTYPE_DIR}"

# ============================================================================
# MERGE MODE: consolidate all genotypes (run after array completes)
# Usage: sbatch --export=MERGE=1 --array=0 run_hla_typing.sbatch
# ============================================================================
if [[ "${MERGE:-0}" == "1" ]]; then
    log_step "hla_typing" "Running arcasHLA merge to consolidate genotypes"

    MERGED_DIR="${OUTDIR}/merged"
    mkdir -p "${MERGED_DIR}"

    # Merge all genotype results
    arcasHLA merge \
        --outdir "${MERGED_DIR}" \
        "${GENOTYPE_DIR}"/*.genotype.json

    log_step "hla_typing" "Merged HLA genotypes written to: ${MERGED_DIR}"
    echo "Merge complete - $(date)"
    exit 0
fi

# ============================================================================
# PER-SAMPLE MODE: extract and genotype
# ============================================================================

# --- Determine sample for this array task ---
mapfile -t SAMPLE_IDS < <(get_sample_ids)
SAMPLE_ID="${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}"
BAM_FILE=$(get_bam_for_sample "${SAMPLE_ID}")

log_step "hla_typing" "Starting HLA typing for sample: ${SAMPLE_ID}"

# --- Validate inputs ---
check_file_exists "${BAM_FILE}" "BAM file for ${SAMPLE_ID}"

# --- Step 1: arcasHLA extract (extract HLA reads from BAM) ---
log_step "hla_typing" "[${SAMPLE_ID}] Extracting HLA reads"

SAMPLE_EXTRACT_DIR="${EXTRACT_DIR}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_EXTRACT_DIR}"

arcasHLA extract \
    --outdir "${SAMPLE_EXTRACT_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --verbose \
    "${BAM_FILE}"

log_step "hla_typing" "[${SAMPLE_ID}] HLA read extraction complete"

# --- Step 2: arcasHLA genotype (type class I + II alleles) ---
log_step "hla_typing" "[${SAMPLE_ID}] Genotyping HLA alleles (class I + II)"

# Find the extracted FASTQ files
EXTRACTED_FQ1="${SAMPLE_EXTRACT_DIR}/${SAMPLE_ID}.extracted.1.fq.gz"
EXTRACTED_FQ2="${SAMPLE_EXTRACT_DIR}/${SAMPLE_ID}.extracted.2.fq.gz"

# arcasHLA genotype with class I and II genes
arcasHLA genotype \
    --outdir "${GENOTYPE_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --genes "${HLA_GENES}" \
    --verbose \
    "${EXTRACTED_FQ1}" "${EXTRACTED_FQ2}"

log_step "hla_typing" "[${SAMPLE_ID}] HLA genotyping complete"
log_step "hla_typing" "[${SAMPLE_ID}] Results: ${GENOTYPE_DIR}/${SAMPLE_ID}.genotype.json"

echo "Done: ${SAMPLE_ID} - $(date)"
