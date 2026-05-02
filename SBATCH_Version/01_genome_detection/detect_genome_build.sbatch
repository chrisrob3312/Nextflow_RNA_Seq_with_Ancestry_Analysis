#!/bin/bash
# ============================================================================
# DETECT GENOME BUILD FROM BAM HEADERS
# ============================================================================
# Determines whether each sample is aligned to hg38 or hg19 by inspecting
# BAM headers. If hg19 samples are found and a chain file is provided,
# performs liftover with CrossMap.
#
# Submit: sbatch detect_genome_build.sbatch
# ============================================================================

#SBATCH --job-name=detect_genome_build
#SBATCH --output=../logs/genome_detection/detect_build_%j.out
#SBATCH --error=../logs/genome_detection/detect_build_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
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
# module load crossmap/0.6.5
# module load python/3.10
# conda activate rnaseq_pipeline

# ============================================================================
# SETUP
# ============================================================================

OUTDIR="${RESULTS_DIR}/genome_detection"
mkdir -p "${OUTDIR}"

log_step "genome_detection" "Starting genome build detection"

# Expected chr1 contig lengths
HG38_CHR1_LENGTH=248956422
HG19_CHR1_LENGTH=249250621

# Summary file
SUMMARY_TSV="${OUTDIR}/genome_build_summary.tsv"
echo -e "sample_id\tdetected_build\tchr_prefix\tchr1_length\tneeds_liftover" > "${SUMMARY_TSV}"

# Track hg19 samples for potential liftover
declare -a HG19_SAMPLES=()

# ============================================================================
# DETECT BUILD FOR EACH SAMPLE
# ============================================================================

log_step "genome_detection" "Processing samples from samplesheet: ${SAMPLESHEET}"

while IFS=',' read -r sample_id bam_path bai_path rest; do
    [[ -z "${sample_id}" ]] && continue

    echo "Processing: ${sample_id}"

    # Check BAM exists
    if [[ ! -f "${bam_path}" ]]; then
        echo "  WARNING: BAM not found, skipping: ${bam_path}"
        echo -e "${sample_id}\tUNKNOWN\tUNKNOWN\t0\tfalse" >> "${SUMMARY_TSV}"
        continue
    fi

    # Extract header
    HEADER=$(samtools view -H "${bam_path}" 2>/dev/null)

    # Detect chromosome naming convention (chr prefix vs numeric)
    CHR_PREFIX="unknown"
    if echo "${HEADER}" | grep -q "SN:chr1"; then
        CHR_PREFIX="chr"
    elif echo "${HEADER}" | grep -q "SN:1"; then
        CHR_PREFIX="numeric"
    fi

    # Get chr1 contig length
    CHR1_LENGTH=0
    if [[ "${CHR_PREFIX}" == "chr" ]]; then
        CHR1_LENGTH=$(echo "${HEADER}" | grep -P "SN:chr1\t" | sed 's/.*LN://' | cut -f1)
    elif [[ "${CHR_PREFIX}" == "numeric" ]]; then
        CHR1_LENGTH=$(echo "${HEADER}" | grep -P "SN:1\t" | sed 's/.*LN://' | cut -f1)
    fi

    # Determine build
    DETECTED_BUILD="unknown"
    NEEDS_LIFTOVER="false"

    if [[ "${CHR1_LENGTH}" -eq "${HG38_CHR1_LENGTH}" ]]; then
        DETECTED_BUILD="hg38"
    elif [[ "${CHR1_LENGTH}" -eq "${HG19_CHR1_LENGTH}" ]]; then
        DETECTED_BUILD="hg19"
        NEEDS_LIFTOVER="true"
        HG19_SAMPLES+=("${sample_id}")
    else
        echo "  WARNING: Unrecognized chr1 length (${CHR1_LENGTH}) for ${sample_id}"
        DETECTED_BUILD="unknown"
    fi

    # Write per-sample result
    SAMPLE_RESULT="${OUTDIR}/${sample_id}.build.txt"
    cat > "${SAMPLE_RESULT}" <<EOF
sample_id=${sample_id}
detected_build=${DETECTED_BUILD}
chr_prefix=${CHR_PREFIX}
chr1_length=${CHR1_LENGTH}
needs_liftover=${NEEDS_LIFTOVER}
bam_path=${bam_path}
detection_date=$(date '+%Y-%m-%d %H:%M:%S')
EOF

    # Append to summary
    echo -e "${sample_id}\t${DETECTED_BUILD}\t${CHR_PREFIX}\t${CHR1_LENGTH}\t${NEEDS_LIFTOVER}" >> "${SUMMARY_TSV}"

    echo "  Build: ${DETECTED_BUILD} | Chr prefix: ${CHR_PREFIX} | Chr1 length: ${CHR1_LENGTH}"

done < <(tail -n +2 "${SAMPLESHEET}")

# ============================================================================
# CROSSMAP LIFTOVER FOR HG19 SAMPLES
# ============================================================================

if [[ ${#HG19_SAMPLES[@]} -gt 0 ]]; then
    echo ""
    echo "============================================================================"
    echo "  Found ${#HG19_SAMPLES[@]} hg19 sample(s) requiring liftover"
    echo "============================================================================"

    if [[ -f "${CHAIN_FILE}" ]]; then
        log_step "genome_detection" "Running CrossMap liftover for ${#HG19_SAMPLES[@]} hg19 samples"

        LIFTOVER_DIR="${OUTDIR}/liftover"
        mkdir -p "${LIFTOVER_DIR}"

        for sample_id in "${HG19_SAMPLES[@]}"; do
            BAM_PATH=$(get_bam_for_sample "${sample_id}")
            OUTPUT_BAM="${LIFTOVER_DIR}/${sample_id}.hg38_lifted.bam"

            echo "  Lifting over: ${sample_id}"
            echo "    Input:  ${BAM_PATH}"
            echo "    Output: ${OUTPUT_BAM}"

            # Run CrossMap
            CrossMap.py bam \
                "${CHAIN_FILE}" \
                "${BAM_PATH}" \
                "${OUTPUT_BAM}" \
                -a "${FASTA}" 2>&1 | tee "${LIFTOVER_DIR}/${sample_id}.crossmap.log"

            # Sort the lifted BAM
            SORTED_BAM="${LIFTOVER_DIR}/${sample_id}.hg38_lifted.sorted.bam"
            samtools sort \
                -@ "${SLURM_CPUS_PER_TASK:-2}" \
                -o "${SORTED_BAM}" \
                "${OUTPUT_BAM}"

            # Index the sorted BAM
            samtools index "${SORTED_BAM}"

            # Remove unsorted intermediate
            rm -f "${OUTPUT_BAM}"

            # Update the per-sample build file with liftover info
            SAMPLE_RESULT="${OUTDIR}/${sample_id}.build.txt"
            echo "liftover_status=completed" >> "${SAMPLE_RESULT}"
            echo "lifted_bam=${SORTED_BAM}" >> "${SAMPLE_RESULT}"

            echo "    Liftover complete."
        done
    else
        echo "  WARNING: Chain file not found: ${CHAIN_FILE}"
        echo "  WARNING: Cannot perform liftover. hg19 samples will NOT be lifted."
        echo "  WARNING: Downstream analysis may fail for these samples."

        for sample_id in "${HG19_SAMPLES[@]}"; do
            SAMPLE_RESULT="${OUTDIR}/${sample_id}.build.txt"
            echo "liftover_status=skipped_no_chain_file" >> "${SAMPLE_RESULT}"
        done
    fi
else
    echo ""
    echo "All samples are hg38 (or unknown). No liftover required."
fi

# ============================================================================
# FINAL SUMMARY
# ============================================================================

echo ""
echo "============================================================================"
echo "  GENOME BUILD DETECTION SUMMARY"
echo "============================================================================"

TOTAL=$(tail -n +2 "${SUMMARY_TSV}" | wc -l)
N_HG38=$(tail -n +2 "${SUMMARY_TSV}" | awk -F'\t' '$2=="hg38"' | wc -l)
N_HG19=$(tail -n +2 "${SUMMARY_TSV}" | awk -F'\t' '$2=="hg19"' | wc -l)
N_UNKNOWN=$(tail -n +2 "${SUMMARY_TSV}" | awk -F'\t' '$2=="unknown"' | wc -l)

echo "  Total samples:   ${TOTAL}"
echo "  hg38:            ${N_HG38}"
echo "  hg19:            ${N_HG19}"
echo "  Unknown:         ${N_UNKNOWN}"
echo ""
echo "  Summary file:    ${SUMMARY_TSV}"
echo "  Per-sample dir:  ${OUTDIR}/"
echo "============================================================================"

log_step "genome_detection" "Genome build detection complete (hg38=${N_HG38}, hg19=${N_HG19}, unknown=${N_UNKNOWN})"

if [[ "${N_UNKNOWN}" -gt 0 ]]; then
    echo "WARNING: ${N_UNKNOWN} sample(s) have unknown genome build. Review manually."
    exit 1
fi

exit 0
