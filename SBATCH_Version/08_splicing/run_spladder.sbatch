#!/bin/bash
# ============================================================================
# DIFFERENTIAL SPLICING ANALYSIS WITH SPLADDER
# ============================================================================
# Runs SplAdder to detect alternative splicing events from RNA-seq alignments
# and tests for differential splicing between conditions.
# ============================================================================

#SBATCH --job-name=spladder
#SBATCH --output=logs/spladder_%j.out
#SBATCH --error=logs/spladder_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --account=your_account

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load spladder/3.0.4
# module load python/3.10
# module load samtools/1.18

log_step "spladder" "Starting SplAdder splicing analysis"

# --- Create output directories ---
OUTDIR="${RESULTS_DIR}/splicing/spladder"
mkdir -p "${OUTDIR}"

# --- Build BAM lists for each condition ---
BAMS_ALL=""
BAMS_RELAPSE=""
BAMS_NO_RELAPSE=""

while IFS=',' read -r sample_id bam bai condition rest; do
    [[ "${sample_id}" == "sample_id" ]] && continue  # skip header
    check_file_exists "${bam}" "BAM file for ${sample_id}"

    if [[ -n "${BAMS_ALL}" ]]; then
        BAMS_ALL="${BAMS_ALL},${bam}"
    else
        BAMS_ALL="${bam}"
    fi

    if [[ "${condition}" == "relapse" ]]; then
        if [[ -n "${BAMS_RELAPSE}" ]]; then
            BAMS_RELAPSE="${BAMS_RELAPSE},${bam}"
        else
            BAMS_RELAPSE="${bam}"
        fi
    elif [[ "${condition}" == "no_relapse" ]]; then
        if [[ -n "${BAMS_NO_RELAPSE}" ]]; then
            BAMS_NO_RELAPSE="${BAMS_NO_RELAPSE},${bam}"
        else
            BAMS_NO_RELAPSE="${bam}"
        fi
    fi
done < "${SAMPLESHEET}"

# --- Validate inputs ---
check_file_exists "${GTF}" "GTF annotation"

if [[ -z "${BAMS_RELAPSE}" ]]; then
    echo "ERROR: No BAM files found for relapse group" >&2
    exit 1
fi
if [[ -z "${BAMS_NO_RELAPSE}" ]]; then
    echo "ERROR: No BAM files found for no_relapse group" >&2
    exit 1
fi

log_step "spladder" "Total BAMs: $(echo "${BAMS_ALL}" | tr ',' '\n' | wc -l)"
log_step "spladder" "Relapse BAMs: $(echo "${BAMS_RELAPSE}" | tr ',' '\n' | wc -l)"
log_step "spladder" "No-relapse BAMs: $(echo "${BAMS_NO_RELAPSE}" | tr ',' '\n' | wc -l)"

# --- Step 1: Build splicing graphs ---
log_step "spladder" "Building splicing graphs"

spladder build \
    -b "${BAMS_ALL}" \
    -o "${OUTDIR}" \
    -a "${GTF}" \
    --merge-strat merge_graphs \
    --confidence "${SPLADDER_CONFIDENCE}" \
    --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip

# --- Step 2: Test for differential splicing ---
log_step "spladder" "Testing for differential splicing between conditions"

spladder test \
    -o "${OUTDIR}" \
    --conditionA "${BAMS_RELAPSE}" \
    --conditionB "${BAMS_NO_RELAPSE}" \
    --event_type all

log_step "spladder" "SplAdder analysis complete"
log_step "spladder" "Output directory: ${OUTDIR}"

echo "Done: SplAdder - $(date)"
