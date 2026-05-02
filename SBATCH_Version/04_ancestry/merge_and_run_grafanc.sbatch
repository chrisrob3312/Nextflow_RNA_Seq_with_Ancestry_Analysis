#!/bin/bash
# ============================================================================
# MERGE GRAF SNP VCFs AND RUN GRAF-ANC ANCESTRY INFERENCE
# ============================================================================
# Runs after extract_graf_snps array job completes. Merges per-sample VCFs
# into a multi-sample VCF, runs grafpop, and parses output into TSV tables.
# ============================================================================

#SBATCH --job-name=merge_run_grafanc
#SBATCH --output=logs/merge_run_grafanc_%j.out
#SBATCH --error=logs/merge_run_grafanc_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal
#SBATCH --account=your_account

# --- Source shared configuration ---
source "$(dirname "$0")/../config.sh"

set -euo pipefail

# --- Load required modules ---
# module load bcftools/1.18

# --- Setup ---
log_step "merge_grafanc" "Starting GRAF-anc merge and analysis"

GRAF_SNP_DIR="${RESULTS_DIR}/ancestry/graf_snps"
OUTDIR="${RESULTS_DIR}/ancestry/grafanc"
mkdir -p "${OUTDIR}"

# --- Validate inputs ---
check_file_exists "${GRAFANC_BINARY}" "GRAF-anc binary (grafpop)"

# --- Build list of per-sample VCFs ---
VCF_LIST="${OUTDIR}/graf_vcf_list.txt"
> "${VCF_LIST}"

mapfile -t SAMPLE_IDS < <(get_sample_ids)
for SAMPLE_ID in "${SAMPLE_IDS[@]}"; do
    VCF="${GRAF_SNP_DIR}/${SAMPLE_ID}.graf_snps.vcf.gz"
    check_file_exists "${VCF}" "GRAF VCF for ${SAMPLE_ID}"
    echo "${VCF}" >> "${VCF_LIST}"
done

log_step "merge_grafanc" "Found ${#SAMPLE_IDS[@]} sample VCFs to merge"

# --- Merge all per-sample VCFs into one multi-sample VCF ---
MERGED_VCF="${OUTDIR}/all_samples.graf_snps.merged.vcf.gz"

bcftools merge \
    --file-list "${VCF_LIST}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --missing-to-ref \
    -Oz \
    -o "${MERGED_VCF}"

bcftools index --tbi "${MERGED_VCF}"

log_step "merge_grafanc" "Merged VCF created: ${MERGED_VCF}"

# --- Run GRAF-anc (grafpop) ---
GRAFANC_OUTPUT="${OUTDIR}/grafanc_results.txt"

"${GRAFANC_BINARY}" \
    --vcf "${MERGED_VCF}" \
    --out "${GRAFANC_OUTPUT}" \
    --data "${GRAFANC_DATA}"

log_step "merge_grafanc" "GRAF-anc analysis complete"

# --- Parse GRAF-anc output into structured TSVs ---

# ancestry_proportions.tsv: Sample-level ancestry proportion estimates
PROPORTIONS_FILE="${OUTDIR}/ancestry_proportions.tsv"
echo -e "sample_id\teuropean\tafrican\teast_asian\tamerican\tsouth_asian" > "${PROPORTIONS_FILE}"

# ancestry_categories.tsv: Categorical ancestry assignments
CATEGORIES_FILE="${OUTDIR}/ancestry_categories.tsv"
echo -e "sample_id\tpredicted_ancestry\tconfidence" > "${CATEGORIES_FILE}"

# Parse the GRAF-anc output file
# GRAF-anc output format: Subject, Self-reported, GrafPop, P_e, P_a, P_f
# (column names may vary by version; adjust awk fields as needed)
if [[ -f "${GRAFANC_OUTPUT}" ]]; then
    # Extract proportions (skip header lines starting with #)
    awk -F'\t' '!/^#/ && NR>1 {
        print $1"\t"$4"\t"$6"\t"$5"\t.\t."
    }' "${GRAFANC_OUTPUT}" >> "${PROPORTIONS_FILE}"

    # Extract categorical assignments
    awk -F'\t' '!/^#/ && NR>1 {
        print $1"\t"$3"\t."
    }' "${GRAFANC_OUTPUT}" >> "${CATEGORIES_FILE}"
else
    echo "WARNING: GRAF-anc output not found: ${GRAFANC_OUTPUT}" >&2
fi

log_step "merge_grafanc" "Results written to:"
log_step "merge_grafanc" "  Proportions: ${PROPORTIONS_FILE}"
log_step "merge_grafanc" "  Categories:  ${CATEGORIES_FILE}"

echo "Done - $(date)"
