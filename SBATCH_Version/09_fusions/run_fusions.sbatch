#!/bin/bash
# ============================================================================
# GENE FUSION DETECTION WITH ARRIBA AND FUSIONCATCHER
# ============================================================================
# SLURM array job: runs Arriba and FusionCatcher for each sample to detect
# gene fusions, then merges results from both callers.
# Set --array=0-(N-1) where N is the number of samples.
# After the array completes, submit with SLURM_ARRAY_TASK_ID=merge or run
# the merge step manually.
# ============================================================================

#SBATCH --job-name=fusions
#SBATCH --output=logs/fusions_%A_%a.out
#SBATCH --error=logs/fusions_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --account=your_account
#SBATCH --array=0-N

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load arriba/2.4.0
# module load fusioncatcher/1.33
# module load samtools/1.18
# module load star/2.7.10b
# module load R/4.3.0

# --- Create output directories ---
OUTDIR="${RESULTS_DIR}/fusions"
mkdir -p "${OUTDIR}"

# --- Determine sample for this array task ---
mapfile -t SAMPLE_IDS < <(get_sample_ids)

# --- Check if this is the merge step ---
if [[ "${SLURM_ARRAY_TASK_ID}" == "merge" ]] || \
   [[ "${SLURM_ARRAY_TASK_ID}" -ge "${#SAMPLE_IDS[@]}" ]]; then

    log_step "fusions" "Running merge step - combining results from all callers"

    MERGED_DIR="${OUTDIR}/merged"
    mkdir -p "${MERGED_DIR}"

    # Combine fusion calls from both tools for all samples
    echo -e "sample_id\tfusion_gene\tbreakpoint1\tbreakpoint2\tsplit_reads\tspan_reads\tcaller\tconfidence" \
        > "${MERGED_DIR}/all_fusions_merged.tsv"

    for SAMPLE_ID in "${SAMPLE_IDS[@]}"; do
        ARRIBA_FILE="${OUTDIR}/${SAMPLE_ID}/arriba/fusions.tsv"
        FC_FILE="${OUTDIR}/${SAMPLE_ID}/fusioncatcher/final-list_candidate-fusion-genes.txt"

        # Parse Arriba results
        if [[ -f "${ARRIBA_FILE}" ]]; then
            tail -n +2 "${ARRIBA_FILE}" | awk -v sample="${SAMPLE_ID}" -F'\t' \
                '{print sample"\t"$1"--"$2"\t"$5"\t"$6"\t"$10"\t"$11"\tarriba\t"$15}' \
                >> "${MERGED_DIR}/all_fusions_merged.tsv"
        fi

        # Parse FusionCatcher results
        if [[ -f "${FC_FILE}" ]]; then
            tail -n +2 "${FC_FILE}" | awk -v sample="${SAMPLE_ID}" -F'\t' \
                '{print sample"\t"$1"--"$2"\t"$9"\t"$10"\t"$5"\t"$6"\tfusioncatcher\t"$3}' \
                >> "${MERGED_DIR}/all_fusions_merged.tsv"
        fi
    done

    # Create summary of fusions found by both callers
    awk -F'\t' 'NR>1 {print $1"\t"$2"\t"$7}' "${MERGED_DIR}/all_fusions_merged.tsv" \
        | sort | uniq -c | sort -rn \
        > "${MERGED_DIR}/fusion_summary.txt"

    log_step "fusions" "Merge complete. Output: ${MERGED_DIR}/all_fusions_merged.tsv"
    echo "Done: Fusion merge - $(date)"
    exit 0
fi

# --- Per-sample fusion detection ---
SAMPLE_ID="${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}"
BAM_FILE=$(get_bam_for_sample "${SAMPLE_ID}")

log_step "fusions" "Starting fusion detection for sample: ${SAMPLE_ID}"

# --- Validate inputs ---
check_file_exists "${BAM_FILE}" "BAM file for ${SAMPLE_ID}"
check_file_exists "${GTF}" "GTF annotation"
check_file_exists "${FASTA}" "Reference FASTA"
check_file_exists "${ARRIBA_BLACKLIST}" "Arriba blacklist"
check_file_exists "${ARRIBA_KNOWN_FUSIONS}" "Arriba known fusions"

# --- Create sample output directories ---
SAMPLE_OUTDIR="${OUTDIR}/${SAMPLE_ID}"
ARRIBA_DIR="${SAMPLE_OUTDIR}/arriba"
FC_DIR="${SAMPLE_OUTDIR}/fusioncatcher"
mkdir -p "${ARRIBA_DIR}" "${FC_DIR}"

# ============================================================================
# ARRIBA
# ============================================================================
log_step "fusions" "Running Arriba for ${SAMPLE_ID}"

arriba \
    -x "${BAM_FILE}" \
    -g "${GTF}" \
    -a "${FASTA}" \
    -b "${ARRIBA_BLACKLIST}" \
    -k "${ARRIBA_KNOWN_FUSIONS}" \
    -o "${ARRIBA_DIR}/fusions.tsv" \
    -O "${ARRIBA_DIR}/fusions.discarded.tsv"

# --- Arriba visualization ---
if [[ -s "${ARRIBA_DIR}/fusions.tsv" ]]; then
    log_step "fusions" "Generating Arriba fusion visualizations for ${SAMPLE_ID}"

    draw_fusions.R \
        --fusions "${ARRIBA_DIR}/fusions.tsv" \
        --alignments "${BAM_FILE}" \
        --annotation "${GTF}" \
        --output "${ARRIBA_DIR}/fusions.pdf"
fi

# ============================================================================
# FUSIONCATCHER
# ============================================================================
log_step "fusions" "Running FusionCatcher for ${SAMPLE_ID}"

# Convert BAM to FASTQ for FusionCatcher
READS_DIR="${SAMPLE_OUTDIR}/fastq_tmp"
mkdir -p "${READS_DIR}"

samtools sort -n -@ "${SLURM_CPUS_PER_TASK}" "${BAM_FILE}" \
    | samtools fastq \
        -@ "${SLURM_CPUS_PER_TASK}" \
        -1 "${READS_DIR}/${SAMPLE_ID}_R1.fastq.gz" \
        -2 "${READS_DIR}/${SAMPLE_ID}_R2.fastq.gz" \
        -0 /dev/null \
        -s /dev/null \
        -

# Run FusionCatcher
fusioncatcher \
    -d "${FUSIONCATCHER_DATA}" \
    -i "${READS_DIR}" \
    -o "${FC_DIR}" \
    -p "${SLURM_CPUS_PER_TASK}"

# Clean up temporary FASTQ files
rm -rf "${READS_DIR}"

log_step "fusions" "Fusion detection complete for sample: ${SAMPLE_ID}"
echo "Done: ${SAMPLE_ID} - $(date)"
