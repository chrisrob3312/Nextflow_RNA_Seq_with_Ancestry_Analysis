#!/bin/bash
#SBATCH --job-name=chasmplus
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/chasmplus_%A_%a.out
#SBATCH --error=logs/chasmplus_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-N%10

set -euo pipefail

# ============================================================================
# 05 - CHASMplus Driver Mutation Prediction
# ============================================================================
# Array job: per-sample OpenCRAVAT with CHASMplus annotator for predicting
# driver mutations from filtered VCFs. Identifies variants with significant
# driver potential (CHASMplus p-value < 0.05).
# ============================================================================

# Source shared configuration
source "$(dirname "$0")/../config.sh"

# --- Input parameters ---
VCF_DIR="${VCF_DIR:-${RESULTS_DIR}/variant_calling/filtered_vcfs}"
OUTPUT_DIR="${CHASMPLUS_OUTPUT_DIR:-${RESULTS_DIR}/variant_calling/chasmplus}"
CHASMPLUS_PVALUE="${CHASMPLUS_PVALUE:-0.05}"
CANCER_TYPE="${CANCER_TYPE:-}"

# --- Get sample list ---
mapfile -t SAMPLES < <(get_sample_ids)
NUM_SAMPLES=${#SAMPLES[@]}

SAMPLE_IDX=${SLURM_ARRAY_TASK_ID:-0}
if [[ $SAMPLE_IDX -ge $NUM_SAMPLES ]]; then
    echo "ERROR: Array index ${SAMPLE_IDX} exceeds sample count ${NUM_SAMPLES}"
    exit 1
fi
SAMPLE_ID="${SAMPLES[$SAMPLE_IDX]}"

# --- Validate inputs ---
log_step "CHASMPLUS" "Processing sample: ${SAMPLE_ID} (index ${SAMPLE_IDX})"

# Find VCF for this sample (support both .vcf and .vcf.gz)
SAMPLE_VCF=""
for ext in ".filtered.vcf.gz" ".filtered.vcf" ".vcf.gz" ".vcf"; do
    if [[ -f "${VCF_DIR}/${SAMPLE_ID}${ext}" ]]; then
        SAMPLE_VCF="${VCF_DIR}/${SAMPLE_ID}${ext}"
        break
    fi
done

if [[ -z "${SAMPLE_VCF}" ]]; then
    echo "ERROR: No VCF found for ${SAMPLE_ID} in ${VCF_DIR}"
    exit 1
fi

log_step "CHASMPLUS" "  VCF: ${SAMPLE_VCF}"

# --- Create output directories ---
SAMPLE_OUTDIR="${OUTPUT_DIR}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_OUTDIR}"
mkdir -p "$(dirname "$0")/logs"

# --- Log parameters ---
log_step "CHASMPLUS" "Parameters:"
log_step "CHASMPLUS" "  VCF: ${SAMPLE_VCF}"
log_step "CHASMPLUS" "  Output: ${SAMPLE_OUTDIR}"
log_step "CHASMPLUS" "  P-value threshold: ${CHASMPLUS_PVALUE}"
log_step "CHASMPLUS" "  Cancer type: ${CANCER_TYPE:-auto}"
log_step "CHASMPLUS" "  CPUs: ${SLURM_CPUS_PER_TASK:-4}"

# --- Load software environment ---
# module load opencravat/2.4
# module load python/3.9
# pip install open-cravat
# oc module install chasmplus

# ============================================================================
# 1. Convert VCF to OpenCRAVAT input format if needed
# ============================================================================
log_step "CHASMPLUS" "  Preparing input..."

OC_INPUT="${SAMPLE_OUTDIR}/${SAMPLE_ID}_oc_input.vcf"

if [[ "${SAMPLE_VCF}" == *.gz ]]; then
    # Decompress VCF for OpenCRAVAT
    zcat "${SAMPLE_VCF}" > "${OC_INPUT}"
else
    cp "${SAMPLE_VCF}" "${OC_INPUT}"
fi

# Count variants
N_VARIANTS=$(grep -cv "^#" "${OC_INPUT}" || echo 0)
log_step "CHASMPLUS" "  Input variants: ${N_VARIANTS}"

if [[ "${N_VARIANTS}" -eq 0 ]]; then
    log_step "CHASMPLUS" "  WARNING: No variants in VCF for ${SAMPLE_ID}, skipping"
    echo -e "sample\tchrom\tpos\tref\talt\tgene\taachange\tchasmplus_score\tchasmplus_pvalue\tdriver_call" \
        > "${SAMPLE_OUTDIR}/${SAMPLE_ID}_driver_mutations.tsv"
    log_step "CHASMPLUS" "  Created empty output file"
    exit 0
fi

# ============================================================================
# 2. Run OpenCRAVAT with CHASMplus annotator
# ============================================================================
log_step "CHASMPLUS" "  Running OpenCRAVAT with CHASMplus annotator..."

# Build annotator list
ANNOTATORS="chasmplus"
if [[ -n "${CANCER_TYPE}" ]]; then
    # Use cancer-type-specific CHASMplus if available
    ANNOTATORS="chasmplus chasmplus_${CANCER_TYPE}"
fi

oc run "${OC_INPUT}" \
    -t vcf \
    -a ${ANNOTATORS} \
    -l hg38 \
    -d "${SAMPLE_OUTDIR}" \
    --mp "${SLURM_CPUS_PER_TASK:-4}" \
    -n "${SAMPLE_ID}" \
    --skip-variant-deduplication

log_step "CHASMPLUS" "  OpenCRAVAT run complete"

# ============================================================================
# 3. Parse OpenCRAVAT output and extract driver predictions
# ============================================================================
log_step "CHASMPLUS" "  Parsing CHASMplus results..."

OC_TSV="${SAMPLE_OUTDIR}/${SAMPLE_ID}.variant.tsv"
OC_SQLITE="${SAMPLE_OUTDIR}/${SAMPLE_ID}.sqlite"

# Try to extract results from TSV output
DRIVER_FILE="${SAMPLE_OUTDIR}/${SAMPLE_ID}_driver_mutations.tsv"

if [[ -f "${OC_TSV}" ]]; then
    # Parse the variant-level TSV file
    # OpenCRAVAT TSV has comment lines starting with # and a header row
    python3 - "${OC_TSV}" "${DRIVER_FILE}" "${SAMPLE_ID}" "${CHASMPLUS_PVALUE}" <<'PYEOF'
import sys
import csv

input_file = sys.argv[1]
output_file = sys.argv[2]
sample_id = sys.argv[3]
pval_thresh = float(sys.argv[4])

# Read OpenCRAVAT output
header = None
rows = []
with open(input_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        if header is None:
            header = line.strip().split('\t')
            continue
        rows.append(line.strip().split('\t'))

if not header or not rows:
    # Write empty output
    with open(output_file, 'w') as out:
        out.write("sample\tchrom\tpos\tref\talt\tgene\taachange\t"
                  "chasmplus_score\tchasmplus_pvalue\tdriver_call\n")
    sys.exit(0)

# Find relevant column indices
col_map = {h.lower(): i for i, h in enumerate(header)}

# Try multiple possible column names for CHASMplus fields
score_cols = ['chasmplus__score', 'chasmplus_score', 'chasmplus.score']
pval_cols = ['chasmplus__pvalue', 'chasmplus_pvalue', 'chasmplus.pvalue',
             'chasmplus__all_pvalue', 'chasmplus_all_pvalue']
chrom_cols = ['chrom', 'chromosome', 'chr']
pos_cols = ['pos', 'position']
ref_cols = ['ref_base', 'ref', 'reference']
alt_cols = ['alt_base', 'alt', 'alternate']
gene_cols = ['hugo', 'gene', 'base__hugo', 'base.hugo']
aa_cols = ['achange', 'aachange', 'base__achange', 'base.achange',
           'protein_change']

def find_col(candidates):
    for c in candidates:
        if c in col_map:
            return col_map[c]
    return None

score_idx = find_col(score_cols)
pval_idx = find_col(pval_cols)
chrom_idx = find_col(chrom_cols)
pos_idx = find_col(pos_cols)
ref_idx = find_col(ref_cols)
alt_idx = find_col(alt_cols)
gene_idx = find_col(gene_cols)
aa_idx = find_col(aa_cols)

def safe_get(row, idx, default=''):
    if idx is None or idx >= len(row):
        return default
    return row[idx]

# Write driver mutations
with open(output_file, 'w') as out:
    out.write("sample\tchrom\tpos\tref\talt\tgene\taachange\t"
              "chasmplus_score\tchasmplus_pvalue\tdriver_call\n")

    for row in rows:
        score = safe_get(row, score_idx, '')
        pval = safe_get(row, pval_idx, '')

        try:
            score_val = float(score) if score else None
            pval_val = float(pval) if pval else None
        except ValueError:
            score_val = None
            pval_val = None

        if pval_val is not None and pval_val < pval_thresh:
            driver = "DRIVER"
        elif pval_val is not None:
            driver = "PASSENGER"
        else:
            driver = "NA"

        # Only output variants with CHASMplus scores
        if score_val is not None:
            out.write(f"{sample_id}\t"
                      f"{safe_get(row, chrom_idx)}\t"
                      f"{safe_get(row, pos_idx)}\t"
                      f"{safe_get(row, ref_idx)}\t"
                      f"{safe_get(row, alt_idx)}\t"
                      f"{safe_get(row, gene_idx)}\t"
                      f"{safe_get(row, aa_idx)}\t"
                      f"{score}\t{pval}\t{driver}\n")

print(f"Parsed CHASMplus results for {sample_id}")
PYEOF

elif [[ -f "${OC_SQLITE}" ]]; then
    # Fall back to SQLite output if TSV not available
    log_step "CHASMPLUS" "  Parsing from SQLite database..."

    python3 - "${OC_SQLITE}" "${DRIVER_FILE}" "${SAMPLE_ID}" "${CHASMPLUS_PVALUE}" <<'PYEOF'
import sys
import sqlite3

db_file = sys.argv[1]
output_file = sys.argv[2]
sample_id = sys.argv[3]
pval_thresh = float(sys.argv[4])

conn = sqlite3.connect(db_file)
cursor = conn.cursor()

# Get table names to find CHASMplus results
cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
tables = [row[0] for row in cursor.fetchall()]

with open(output_file, 'w') as out:
    out.write("sample\tchrom\tpos\tref\talt\tgene\taachange\t"
              "chasmplus_score\tchasmplus_pvalue\tdriver_call\n")

    try:
        # Query variant table with CHASMplus annotations
        query = """
            SELECT v.base__chrom, v.base__pos, v.base__ref_base,
                   v.base__alt_base, v.base__hugo, v.base__achange,
                   c.score, c.pvalue
            FROM variant v
            LEFT JOIN variant_chasmplus c ON v.base__uid = c.base__uid
            WHERE c.score IS NOT NULL
        """
        cursor.execute(query)

        for row in cursor.fetchall():
            chrom, pos, ref, alt, gene, aa, score, pval = row
            if pval is not None and pval < pval_thresh:
                driver = "DRIVER"
            elif pval is not None:
                driver = "PASSENGER"
            else:
                driver = "NA"
            out.write(f"{sample_id}\t{chrom}\t{pos}\t{ref}\t{alt}\t"
                      f"{gene}\t{aa}\t{score}\t{pval}\t{driver}\n")
    except Exception as e:
        print(f"WARNING: SQLite query failed: {e}")

conn.close()
print(f"Parsed CHASMplus results from SQLite for {sample_id}")
PYEOF

else
    log_step "CHASMPLUS" "  WARNING: No OpenCRAVAT output found for ${SAMPLE_ID}"
    echo -e "sample\tchrom\tpos\tref\talt\tgene\taachange\tchasmplus_score\tchasmplus_pvalue\tdriver_call" \
        > "${DRIVER_FILE}"
fi

# ============================================================================
# 4. Summary statistics
# ============================================================================
if [[ -f "${DRIVER_FILE}" ]]; then
    N_SCORED=$(tail -n +2 "${DRIVER_FILE}" | wc -l)
    N_DRIVERS=$(tail -n +2 "${DRIVER_FILE}" | awk -F'\t' '$10=="DRIVER"' | wc -l)
    N_PASSENGERS=$(tail -n +2 "${DRIVER_FILE}" | awk -F'\t' '$10=="PASSENGER"' | wc -l)

    log_step "CHASMPLUS" "  Results for ${SAMPLE_ID}:"
    log_step "CHASMPLUS" "    Scored variants: ${N_SCORED}"
    log_step "CHASMPLUS" "    Predicted drivers (p < ${CHASMPLUS_PVALUE}): ${N_DRIVERS}"
    log_step "CHASMPLUS" "    Predicted passengers: ${N_PASSENGERS}"

    # Extract driver gene list
    if [[ ${N_DRIVERS} -gt 0 ]]; then
        DRIVER_GENES=$(tail -n +2 "${DRIVER_FILE}" | \
            awk -F'\t' '$10=="DRIVER" {print $6}' | sort -u | tr '\n' ',' | sed 's/,$//')
        log_step "CHASMPLUS" "    Driver genes: ${DRIVER_GENES}"
    fi
fi

# --- Clean up ---
rm -f "${OC_INPUT}"

# --- Summary ---
echo ""
echo "============================================"
echo "CHASMplus Summary: ${SAMPLE_ID}"
echo "============================================"
echo "  Input VCF: ${SAMPLE_VCF}"
echo "  Input variants: ${N_VARIANTS}"
echo "  Scored variants: ${N_SCORED:-0}"
echo "  Predicted drivers: ${N_DRIVERS:-0}"
echo "  Output: ${DRIVER_FILE}"
echo "============================================"

log_step "CHASMPLUS" "Sample ${SAMPLE_ID} CHASMplus analysis complete"
