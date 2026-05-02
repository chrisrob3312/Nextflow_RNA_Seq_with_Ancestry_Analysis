#!/bin/bash
# ============================================================================
# DIFFERENTIAL SPLICING ANALYSIS WITH rMATS TURBO
# ============================================================================
# Runs rMATS for differential alternative splicing between two groups
# (relapse vs no_relapse) defined in the sample metadata.
# ============================================================================

#SBATCH --job-name=rmats
#SBATCH --output=logs/rmats_%j.out
#SBATCH --error=logs/rmats_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --account=your_account

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load rmats/4.3.0
# module load python/3.10
# module load star/2.7.10b

log_step "rmats" "Starting rMATS differential splicing analysis"

# --- Create output directories ---
OUTDIR="${RESULTS_DIR}/splicing/rmats"
TMPDIR_RMATS="${TEMP_DIR}/rmats_tmp"
mkdir -p "${OUTDIR}" "${TMPDIR_RMATS}"

# --- Build BAM list files from metadata ---
# Samplesheet expected format: sample_id,bam,bai,condition,...
# Field 4 is assumed to be the condition column (relapse / no_relapse)
GROUP1_FILE="${OUTDIR}/group1.txt"
GROUP2_FILE="${OUTDIR}/group2.txt"

> "${GROUP1_FILE}"
> "${GROUP2_FILE}"

while IFS=',' read -r sample_id bam bai condition rest; do
    [[ "${sample_id}" == "sample_id" ]] && continue  # skip header
    check_file_exists "${bam}" "BAM file for ${sample_id}"
    if [[ "${condition}" == "relapse" ]]; then
        echo "${bam}" >> "${GROUP1_FILE}"
    elif [[ "${condition}" == "no_relapse" ]]; then
        echo "${bam}" >> "${GROUP2_FILE}"
    fi
done < "${SAMPLESHEET}"

# Verify we have samples in both groups
if [[ ! -s "${GROUP1_FILE}" ]]; then
    echo "ERROR: No samples found in relapse group" >&2
    exit 1
fi
if [[ ! -s "${GROUP2_FILE}" ]]; then
    echo "ERROR: No samples found in no_relapse group" >&2
    exit 1
fi

log_step "rmats" "Group 1 (relapse): $(wc -l < "${GROUP1_FILE}") samples"
log_step "rmats" "Group 2 (no_relapse): $(wc -l < "${GROUP2_FILE}") samples"

# --- Validate reference ---
check_file_exists "${GTF}" "GTF annotation"

# --- Run rMATS turbo ---
rmats.py \
    --b1 "${GROUP1_FILE}" \
    --b2 "${GROUP2_FILE}" \
    --gtf "${GTF}" \
    -t paired \
    --readLength "${RMATS_READ_LENGTH}" \
    --nthread "${SLURM_CPUS_PER_TASK}" \
    --od "${OUTDIR}" \
    --tmp "${TMPDIR_RMATS}" \
    --novelSS

# --- Clean up temp directory ---
rm -rf "${TMPDIR_RMATS}"

log_step "rmats" "rMATS analysis complete"
log_step "rmats" "Output directory: ${OUTDIR}"

echo "Done: rMATS - $(date)"
