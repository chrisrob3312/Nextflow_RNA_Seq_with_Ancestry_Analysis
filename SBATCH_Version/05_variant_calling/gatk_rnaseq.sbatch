#!/bin/bash
# ============================================================================
# GATK RNA-SEQ VARIANT CALLING PIPELINE
# ============================================================================
# SLURM array job: per-sample GATK Best Practices for RNA-seq variant calling.
# Steps: SplitNCigarReads -> BQSR (optional) -> HaplotypeCaller -> Filter -> TMB
# ============================================================================

#SBATCH --job-name=gatk_rnaseq
#SBATCH --output=logs/gatk_rnaseq_%A_%a.out
#SBATCH --error=logs/gatk_rnaseq_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --account=your_account
#SBATCH --array=0-N

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load gatk/4.4.0
# module load samtools/1.18

# --- Determine sample for this array task ---
mapfile -t SAMPLE_IDS < <(get_sample_ids)
SAMPLE_ID="${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}"
BAM_FILE=$(get_bam_for_sample "${SAMPLE_ID}")

log_step "gatk_rnaseq" "Starting GATK RNA-seq variant calling for sample: ${SAMPLE_ID}"

# --- Validate inputs ---
check_file_exists "${BAM_FILE}" "BAM file for ${SAMPLE_ID}"
check_file_exists "${FASTA}" "Reference FASTA"
check_file_exists "${FASTA_FAI}" "Reference FASTA index"

# --- Create output directories ---
OUTDIR="${RESULTS_DIR}/variant_calling/${SAMPLE_ID}"
TMPDIR_SAMPLE="${TEMP_DIR}/gatk_${SAMPLE_ID}"
mkdir -p "${OUTDIR}" "${TMPDIR_SAMPLE}"

# ============================================================================
# Step 1: SplitNCigarReads
# ============================================================================
log_step "gatk_rnaseq" "[${SAMPLE_ID}] Running SplitNCigarReads"

SPLIT_BAM="${TMPDIR_SAMPLE}/${SAMPLE_ID}.split.bam"

gatk SplitNCigarReads \
    -R "${FASTA}" \
    -I "${BAM_FILE}" \
    -O "${SPLIT_BAM}" \
    --tmp-dir "${TMPDIR_SAMPLE}" \
    --java-options "-Xmx28g -Djava.io.tmpdir=${TMPDIR_SAMPLE}"

# ============================================================================
# Step 2: Base Quality Score Recalibration (if known_snps provided)
# ============================================================================
RECAL_BAM="${SPLIT_BAM}"

if [[ -n "${KNOWN_SNPS:-}" && -f "${KNOWN_SNPS}" ]]; then
    log_step "gatk_rnaseq" "[${SAMPLE_ID}] Running BaseRecalibrator"

    RECAL_TABLE="${TMPDIR_SAMPLE}/${SAMPLE_ID}.recal_data.table"

    gatk BaseRecalibrator \
        -R "${FASTA}" \
        -I "${SPLIT_BAM}" \
        --known-sites "${KNOWN_SNPS}" \
        -O "${RECAL_TABLE}" \
        --tmp-dir "${TMPDIR_SAMPLE}" \
        --java-options "-Xmx28g -Djava.io.tmpdir=${TMPDIR_SAMPLE}"

    log_step "gatk_rnaseq" "[${SAMPLE_ID}] Running ApplyBQSR"

    RECAL_BAM="${TMPDIR_SAMPLE}/${SAMPLE_ID}.recal.bam"

    gatk ApplyBQSR \
        -R "${FASTA}" \
        -I "${SPLIT_BAM}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${RECAL_BAM}" \
        --tmp-dir "${TMPDIR_SAMPLE}" \
        --java-options "-Xmx28g -Djava.io.tmpdir=${TMPDIR_SAMPLE}"
else
    log_step "gatk_rnaseq" "[${SAMPLE_ID}] Skipping BQSR (no known_snps provided)"
fi

# ============================================================================
# Step 3: HaplotypeCaller (RNA-seq mode)
# ============================================================================
log_step "gatk_rnaseq" "[${SAMPLE_ID}] Running HaplotypeCaller"

RAW_VCF="${OUTDIR}/${SAMPLE_ID}.raw.vcf.gz"

gatk HaplotypeCaller \
    -R "${FASTA}" \
    -I "${RECAL_BAM}" \
    -O "${RAW_VCF}" \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20.0 \
    --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK}" \
    --tmp-dir "${TMPDIR_SAMPLE}" \
    --java-options "-Xmx28g -Djava.io.tmpdir=${TMPDIR_SAMPLE}"

# ============================================================================
# Step 4: VariantFiltration
# ============================================================================
log_step "gatk_rnaseq" "[${SAMPLE_ID}] Running VariantFiltration"

FILTERED_VCF="${OUTDIR}/${SAMPLE_ID}.filtered.vcf.gz"

gatk VariantFiltration \
    -R "${FASTA}" \
    -V "${RAW_VCF}" \
    -O "${FILTERED_VCF}" \
    --filter-expression "FS > 30.0" --filter-name "FS_filter" \
    --filter-expression "QD < 2.0" --filter-name "QD_filter" \
    --filter-expression "MQ < 40.0" --filter-name "MQ_filter" \
    --tmp-dir "${TMPDIR_SAMPLE}" \
    --java-options "-Xmx28g -Djava.io.tmpdir=${TMPDIR_SAMPLE}"

# ============================================================================
# Step 5: Simple TMB Calculation
# ============================================================================
log_step "gatk_rnaseq" "[${SAMPLE_ID}] Calculating TMB"

TMB_FILE="${OUTDIR}/${SAMPLE_ID}.tmb.txt"

# Count PASS variants (SNVs only for TMB)
PASS_VARIANTS=$(bcftools view -f PASS -v snps "${FILTERED_VCF}" 2>/dev/null | grep -cv "^#" || echo 0)

# Estimate callable megabases from the BAM (regions with sufficient coverage)
# Using a simple approximation: count bases with >= 10x coverage
CALLABLE_MB=$(samtools depth -Q 20 -q 20 "${BAM_FILE}" | \
    awk '$3 >= 10 {count++} END {printf "%.2f", count/1000000}')

# Calculate TMB (mutations per megabase)
if [[ $(echo "${CALLABLE_MB} > 0" | bc -l) -eq 1 ]]; then
    TMB=$(echo "scale=4; ${PASS_VARIANTS} / ${CALLABLE_MB}" | bc -l)
else
    TMB="NA"
    CALLABLE_MB="0"
fi

# Write TMB results
cat > "${TMB_FILE}" <<EOF
sample_id	pass_snvs	callable_mb	tmb_per_mb
${SAMPLE_ID}	${PASS_VARIANTS}	${CALLABLE_MB}	${TMB}
EOF

log_step "gatk_rnaseq" "[${SAMPLE_ID}] TMB: ${TMB} mutations/Mb (${PASS_VARIANTS} PASS SNVs / ${CALLABLE_MB} Mb)"

# --- Cleanup temporary files ---
rm -rf "${TMPDIR_SAMPLE}"

log_step "gatk_rnaseq" "Completed GATK RNA-seq variant calling for sample: ${SAMPLE_ID}"
log_step "gatk_rnaseq" "Output: ${FILTERED_VCF}"

echo "Done: ${SAMPLE_ID} - $(date)"
