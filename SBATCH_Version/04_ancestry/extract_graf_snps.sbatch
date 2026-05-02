#!/bin/bash
# ============================================================================
# EXTRACT GRAF-ANC SNPs FROM RNA-SEQ BAMs
# ============================================================================
# SLURM array job: runs bcftools mpileup + call at GRAF-anc 282K SNP positions
# for each sample, producing per-sample VCFs for ancestry inference.
# ============================================================================

#SBATCH --job-name=extract_graf_snps
#SBATCH --output=logs/extract_graf_snps_%A_%a.out
#SBATCH --error=logs/extract_graf_snps_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal
#SBATCH --account=your_account
#SBATCH --array=0-N

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load bcftools/1.18
# module load samtools/1.18

# --- Determine sample for this array task ---
mapfile -t SAMPLE_IDS < <(get_sample_ids)
SAMPLE_ID="${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}"
BAM_FILE=$(get_bam_for_sample "${SAMPLE_ID}")

log_step "extract_graf_snps" "Starting GRAF SNP extraction for sample: ${SAMPLE_ID}"

# --- Validate inputs ---
check_file_exists "${BAM_FILE}" "BAM file for ${SAMPLE_ID}"
check_file_exists "${FASTA}" "Reference FASTA"
check_file_exists "${GRAF_SNP_BED}" "GRAF SNP BED file"

# --- Create output directory ---
OUTDIR="${RESULTS_DIR}/ancestry/graf_snps"
mkdir -p "${OUTDIR}"

# --- Extract GRAF-anc SNP positions ---
# The GRAF_SNP_BED contains the 282K ancestry-informative SNP positions.
# We use bcftools mpileup to generate genotype likelihoods and bcftools call
# to produce genotype calls at these specific positions.

OUTPUT_VCF="${OUTDIR}/${SAMPLE_ID}.graf_snps.vcf.gz"

bcftools mpileup \
    --fasta-ref "${FASTA}" \
    --regions-file "${GRAF_SNP_BED}" \
    --annotate FORMAT/AD,FORMAT/DP \
    --max-depth 10000 \
    --min-MQ 20 \
    --min-BQ 20 \
    --threads "${SLURM_CPUS_PER_TASK}" \
    "${BAM_FILE}" \
| bcftools call \
    -m \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -Oz \
    -o "${OUTPUT_VCF}"

# --- Index the output VCF ---
bcftools index --tbi "${OUTPUT_VCF}"

log_step "extract_graf_snps" "Completed GRAF SNP extraction for sample: ${SAMPLE_ID}"
log_step "extract_graf_snps" "Output: ${OUTPUT_VCF}"

echo "Done: ${SAMPLE_ID} - $(date)"
