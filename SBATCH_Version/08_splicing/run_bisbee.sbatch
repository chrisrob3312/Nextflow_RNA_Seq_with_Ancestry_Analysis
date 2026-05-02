#!/bin/bash
# ============================================================================
# BISBEE: SPLICING ANALYSIS FROM SPLADDER OUTPUT
# ============================================================================
# Runs the Bisbee pipeline (prep, diff, prot, outlier) to perform
# beta-binomial differential splicing, protein effect prediction,
# and per-sample outlier detection using SplAdder HDF5 output.
# Must be run AFTER run_spladder.sbatch completes.
# ============================================================================

#SBATCH --job-name=bisbee
#SBATCH --output=logs/bisbee_%j.out
#SBATCH --error=logs/bisbee_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --account=your_account
#SBATCH --dependency=afterok:SPLADDER_JOB_ID

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load python/3.10
# module load R/4.3.0

log_step "bisbee" "Starting Bisbee splicing analysis"

# --- Define paths ---
BIN_DIR="${PROJECT_DIR}/bin"
SPLADDER_DIR="${RESULTS_DIR}/splicing/spladder"
OUTDIR="${RESULTS_DIR}/splicing/bisbee"
mkdir -p "${OUTDIR}"

# --- Validate SplAdder output exists ---
if [[ ! -d "${SPLADDER_DIR}" ]]; then
    echo "ERROR: SplAdder output directory not found: ${SPLADDER_DIR}" >&2
    echo "Ensure run_spladder.sbatch has completed successfully." >&2
    exit 1
fi

# --- Build groups file for bisbee ---
GROUPS_FILE="${OUTDIR}/groups_file.txt"
> "${GROUPS_FILE}"

while IFS=',' read -r sample_id bam bai condition rest; do
    [[ "${sample_id}" == "sample_id" ]] && continue  # skip header
    if [[ "${condition}" == "relapse" || "${condition}" == "no_relapse" ]]; then
        echo "${sample_id},${condition}" >> "${GROUPS_FILE}"
    fi
done < "${SAMPLESHEET}"

# --- Step 1: Prep - extract junction counts from SplAdder HDF5 output ---
log_step "bisbee" "Running bisbee_prep.py - extracting junction counts"

python "${BIN_DIR}/bisbee_prep.py" \
    --spladder_dir "${SPLADDER_DIR}" \
    --output_dir "${OUTDIR}/prep"

# --- Step 2: Diff - beta-binomial differential splicing ---
log_step "bisbee" "Running bisbee_diff.R - differential splicing analysis"

Rscript "${BIN_DIR}/bisbee_diff.R" \
    --counts_dir "${OUTDIR}/prep" \
    --groups_file "${GROUPS_FILE}" \
    --output_dir "${OUTDIR}/diff" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --fdr "${BISBEE_OUTLIER_FDR}"

# --- Step 3: Prot - protein effect prediction ---
log_step "bisbee" "Running bisbee_prot.py - protein effect prediction"

python "${BIN_DIR}/bisbee_prot.py" \
    --events_dir "${OUTDIR}/diff" \
    --gtf "${GTF}" \
    --fasta "${FASTA}" \
    --output_dir "${OUTDIR}/prot"

# --- Step 4: Outlier - per-sample outlier detection ---
log_step "bisbee" "Running bisbee_outlier.R - per-sample outlier detection"

Rscript "${BIN_DIR}/bisbee_outlier.R" \
    --counts_dir "${OUTDIR}/prep" \
    --groups_file "${GROUPS_FILE}" \
    --output_dir "${OUTDIR}/outlier" \
    --fdr "${BISBEE_OUTLIER_FDR}" \
    --threads "${SLURM_CPUS_PER_TASK}"

log_step "bisbee" "Bisbee analysis complete"
log_step "bisbee" "Output directory: ${OUTDIR}"

echo "Done: Bisbee - $(date)"
