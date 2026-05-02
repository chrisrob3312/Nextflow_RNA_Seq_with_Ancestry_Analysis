#!/bin/bash
# ============================================================================
# MASTER PIPELINE SUBMISSION SCRIPT
# ============================================================================
# Submits all cancer RNA-seq pipeline steps with proper SLURM dependency
# chains so that each step runs only after its prerequisites complete.
#
# Usage:
#   ./submit_all.sh
#
# The script will:
#   1. Source config.sh for shared variables
#   2. Validate that the samplesheet and key references exist
#   3. Submit all sbatch jobs with --dependency=afterok chaining
#   4. Print a summary of submitted job IDs
#
# Dependency structure:
#   setup -> genome_detection -> qc (parallel with counting)
#   counting -> ancestry, variant_calling, hla_typing (parallel)
#   ancestry + counting -> DE -> wgcna, pathway (parallel)
#   counting -> splicing, fusions, cnv (parallel)
#   DE + hla + variants -> neoantigen
#   counting -> immune, tcr (parallel)
#   DE + ancestry -> sensitivity, pharmacogenomics, subtyping (parallel)
#   all above -> validation -> visualization
# ============================================================================

set -euo pipefail

# ============================================================================
# SOURCE CONFIG
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "============================================================================"
echo "  CANCER RNA-SEQ PIPELINE - MASTER JOB SUBMISSION"
echo "============================================================================"
echo ""
echo "  Project: ${PROJECT_DIR}"
echo "  Results: ${RESULTS_DIR}"
echo ""

# ============================================================================
# PRE-FLIGHT VALIDATION
# ============================================================================

echo "[PRE-FLIGHT] Validating inputs..."
ERRORS=0

# Check samplesheet
if [[ ! -f "${SAMPLESHEET}" ]]; then
    echo "  ERROR: Samplesheet not found: ${SAMPLESHEET}"
    ERRORS=$((ERRORS + 1))
else
    echo "  OK: Samplesheet found ($(tail -n +2 "${SAMPLESHEET}" | grep -c . || echo 0) samples)"
fi

# Check key reference files
for ref_file in "${FASTA}" "${GTF}" "${KNOWN_SNPS}"; do
    if [[ ! -f "${ref_file}" ]]; then
        echo "  ERROR: Reference file not found: ${ref_file}"
        ERRORS=$((ERRORS + 1))
    fi
done

if [[ "${ERRORS}" -gt 0 ]]; then
    echo ""
    echo "  ABORTING: ${ERRORS} error(s) found. Fix issues above before submitting."
    exit 1
fi

echo "  All pre-flight checks passed."
echo ""

# ============================================================================
# DETERMINE NUMBER OF SAMPLES (for array jobs)
# ============================================================================

N_SAMPLES=$(tail -n +2 "${SAMPLESHEET}" | grep -c . || echo 0)
if [[ "${N_SAMPLES}" -eq 0 ]]; then
    echo "ERROR: No samples found in samplesheet."
    exit 1
fi
ARRAY_MAX=$((N_SAMPLES - 1))
echo "  Samples: ${N_SAMPLES} (array indices 0-${ARRAY_MAX})"
echo ""

# ============================================================================
# HELPER FUNCTION
# ============================================================================

# Submit a job and capture its ID; optionally add --dependency
submit_job() {
    local script="$1"
    local dep_flag="${2:-}"
    local extra_args="${3:-}"

    if [[ ! -f "${script}" ]]; then
        echo "  WARNING: Script not found, skipping: ${script}"
        echo "SKIP"
        return
    fi

    local cmd="sbatch --parsable"
    if [[ -n "${dep_flag}" ]]; then
        cmd="${cmd} --dependency=${dep_flag}"
    fi
    if [[ -n "${extra_args}" ]]; then
        cmd="${cmd} ${extra_args}"
    fi
    cmd="${cmd} ${script}"

    local job_id
    job_id=$(eval "${cmd}")
    echo "${job_id}"
}

# ============================================================================
# STEP 0: SETUP (run locally, not via SLURM)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 0] Running setup (local)..."
bash "${SCRIPT_DIR}/00_setup/setup_directories.sh"
echo ""

# ============================================================================
# STEP 1: GENOME DETECTION
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 1] Submitting genome detection..."
JOB_GENOME=$(submit_job "${SCRIPT_DIR}/01_genome_detection/detect_genome_build.sbatch")
echo "  Job ID: ${JOB_GENOME}"

# ============================================================================
# STEP 2: QC (array job, depends on genome detection)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 2] Submitting QC (array)..."
if [[ "${JOB_GENOME}" != "SKIP" ]]; then
    JOB_QC=$(submit_job "${SCRIPT_DIR}/02_qc/run_qc.sbatch" "afterok:${JOB_GENOME}" "--array=0-${ARRAY_MAX}")
else
    JOB_QC=$(submit_job "${SCRIPT_DIR}/02_qc/run_qc.sbatch" "" "--array=0-${ARRAY_MAX}")
fi
echo "  Job ID: ${JOB_QC}"

# ============================================================================
# STEP 3: COUNTING (array + merge, depends on genome detection)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 3] Submitting counting..."

# Per-sample featureCounts (array)
if [[ "${JOB_GENOME}" != "SKIP" ]]; then
    JOB_COUNT_ARRAY=$(submit_job "${SCRIPT_DIR}/03_counting/run_featurecounts.sbatch" "afterok:${JOB_GENOME}" "--array=0-${ARRAY_MAX}")
else
    JOB_COUNT_ARRAY=$(submit_job "${SCRIPT_DIR}/03_counting/run_featurecounts.sbatch" "" "--array=0-${ARRAY_MAX}")
fi
echo "  Count array Job ID: ${JOB_COUNT_ARRAY}"

# Merge + normalize (depends on all array tasks)
if [[ "${JOB_COUNT_ARRAY}" != "SKIP" ]]; then
    JOB_COUNT_MERGE=$(submit_job "${SCRIPT_DIR}/03_counting/run_featurecounts.sbatch" "afterok:${JOB_COUNT_ARRAY}" "--export=ALL,MODE=merge")
else
    JOB_COUNT_MERGE="SKIP"
fi
echo "  Count merge Job ID: ${JOB_COUNT_MERGE}"

# Use the merge job as the counting dependency for downstream
JOB_COUNTING="${JOB_COUNT_MERGE}"

# ============================================================================
# STEPS 4-6: ANCESTRY, VARIANT CALLING, HLA TYPING (parallel, after counting)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 4] Submitting ancestry analysis..."

# Extract GRAF SNPs (array, depends on counting)
DEP_COUNTING=""
if [[ "${JOB_COUNTING}" != "SKIP" ]]; then
    DEP_COUNTING="afterok:${JOB_COUNTING}"
fi

JOB_GRAF_EXTRACT=$(submit_job "${SCRIPT_DIR}/04_ancestry/extract_graf_snps.sbatch" "${DEP_COUNTING}" "--array=0-${ARRAY_MAX}")
echo "  GRAF extract Job ID: ${JOB_GRAF_EXTRACT}"

# Merge + run GRAFpop (depends on all extract tasks)
if [[ "${JOB_GRAF_EXTRACT}" != "SKIP" ]]; then
    JOB_ANCESTRY=$(submit_job "${SCRIPT_DIR}/04_ancestry/merge_and_run_grafanc.sbatch" "afterok:${JOB_GRAF_EXTRACT}")
else
    JOB_ANCESTRY=$(submit_job "${SCRIPT_DIR}/04_ancestry/merge_and_run_grafanc.sbatch" "${DEP_COUNTING}")
fi
echo "  GRAF merge Job ID: ${JOB_ANCESTRY}"

# Somalier (parallel with GRAF, depends on counting)
JOB_SOMALIER=$(submit_job "${SCRIPT_DIR}/04_ancestry/run_somalier.sbatch" "${DEP_COUNTING}")
echo "  Somalier Job ID: ${JOB_SOMALIER}"

echo ""
echo "[STEP 5] Submitting variant calling..."
JOB_VARIANTS=$(submit_job "${SCRIPT_DIR}/05_variant_calling/gatk_rnaseq.sbatch" "${DEP_COUNTING}" "--array=0-${ARRAY_MAX}")
echo "  Variant calling Job ID: ${JOB_VARIANTS}"

echo ""
echo "[STEP 6] Submitting HLA typing..."
JOB_HLA=$(submit_job "${SCRIPT_DIR}/06_hla_typing/run_hla_typing.sbatch" "${DEP_COUNTING}" "--array=0-${ARRAY_MAX}")
echo "  HLA typing Job ID: ${JOB_HLA}"

# ============================================================================
# STEP 7: DIFFERENTIAL EXPRESSION (depends on counting + ancestry)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 7] Submitting differential expression..."

# Build dependency: needs counting AND ancestry
DE_DEPS=""
if [[ "${JOB_COUNTING}" != "SKIP" && "${JOB_ANCESTRY}" != "SKIP" ]]; then
    DE_DEPS="afterok:${JOB_COUNTING}:${JOB_ANCESTRY}"
elif [[ "${JOB_COUNTING}" != "SKIP" ]]; then
    DE_DEPS="afterok:${JOB_COUNTING}"
elif [[ "${JOB_ANCESTRY}" != "SKIP" ]]; then
    DE_DEPS="afterok:${JOB_ANCESTRY}"
fi

JOB_DE=$(submit_job "${SCRIPT_DIR}/07_differential_expression/run_de.sbatch" "${DE_DEPS}")
echo "  DE Job ID: ${JOB_DE}"

# ============================================================================
# STEPS 8-10: SPLICING, FUSIONS, CNV (parallel, depend on counting)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 8] Submitting splicing analysis..."
JOB_RMATS=$(submit_job "${SCRIPT_DIR}/08_splicing/run_rmats.sbatch" "${DEP_COUNTING}")
JOB_LEAFCUTTER=$(submit_job "${SCRIPT_DIR}/08_splicing/run_leafcutter.sbatch" "${DEP_COUNTING}")
JOB_SPLADDER=$(submit_job "${SCRIPT_DIR}/08_splicing/run_spladder.sbatch" "${DEP_COUNTING}")
JOB_BISBEE=$(submit_job "${SCRIPT_DIR}/08_splicing/run_bisbee.sbatch" "${DEP_COUNTING}")
echo "  rMATS Job ID: ${JOB_RMATS}"
echo "  Leafcutter Job ID: ${JOB_LEAFCUTTER}"
echo "  SpLaDDer Job ID: ${JOB_SPLADDER}"
echo "  BISBEE Job ID: ${JOB_BISBEE}"

echo ""
echo "[STEP 9] Submitting fusion detection..."
JOB_FUSIONS=$(submit_job "${SCRIPT_DIR}/09_fusions/run_fusions.sbatch" "${DEP_COUNTING}" "--array=0-${ARRAY_MAX}")
echo "  Fusions Job ID: ${JOB_FUSIONS}"

echo ""
echo "[STEP 10] Submitting CNV inference..."
JOB_CNV=$(submit_job "${SCRIPT_DIR}/10_cnv_inference/run_infercnv.sbatch" "${DEP_COUNTING}")
echo "  CNV Job ID: ${JOB_CNV}"

# ============================================================================
# STEPS 11-12: WGCNA, PATHWAY (depend on DE)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 11] Submitting WGCNA..."
DEP_DE=""
if [[ "${JOB_DE}" != "SKIP" ]]; then
    DEP_DE="afterok:${JOB_DE}"
fi
JOB_WGCNA=$(submit_job "${SCRIPT_DIR}/11_wgcna/run_wgcna.sbatch" "${DEP_DE}")
echo "  WGCNA Job ID: ${JOB_WGCNA}"

echo ""
echo "[STEP 12] Submitting pathway enrichment..."
JOB_PATHWAY=$(submit_job "${SCRIPT_DIR}/12_pathway_enrichment/run_pathway.sbatch" "${DEP_DE}")
echo "  Pathway Job ID: ${JOB_PATHWAY}"

# ============================================================================
# STEP 13: IMMUNE ANALYSIS (depends on counting)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 13] Submitting immune deconvolution..."
JOB_IMMUNE=$(submit_job "${SCRIPT_DIR}/13_immune_analysis/run_immune.sbatch" "${DEP_COUNTING}")
echo "  Immune Job ID: ${JOB_IMMUNE}"

# ============================================================================
# STEP 14: NEOANTIGEN (depends on DE + HLA + variants)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 14] Submitting neoantigen prediction..."

# Build dependency: needs DE, HLA, and variants
NEO_DEPS_LIST=()
if [[ "${JOB_DE}" != "SKIP" ]]; then NEO_DEPS_LIST+=("${JOB_DE}"); fi
if [[ "${JOB_HLA}" != "SKIP" ]]; then NEO_DEPS_LIST+=("${JOB_HLA}"); fi
if [[ "${JOB_VARIANTS}" != "SKIP" ]]; then NEO_DEPS_LIST+=("${JOB_VARIANTS}"); fi

NEO_DEPS=""
if [[ ${#NEO_DEPS_LIST[@]} -gt 0 ]]; then
    NEO_DEPS="afterok:$(IFS=:; echo "${NEO_DEPS_LIST[*]}")"
fi

JOB_NEOANTIGEN=$(submit_job "${SCRIPT_DIR}/14_neoantigen/run_neoantigen.sbatch" "${NEO_DEPS}")
echo "  Neoantigen Job ID: ${JOB_NEOANTIGEN}"

# ============================================================================
# STEP 15: TCR REPERTOIRE (depends on counting)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 15] Submitting TCR repertoire analysis..."
JOB_TCR=$(submit_job "${SCRIPT_DIR}/15_tcr_repertoire/run_trust4.sbatch" "${DEP_COUNTING}" "--array=0-${ARRAY_MAX}")
echo "  TCR Job ID: ${JOB_TCR}"

# ============================================================================
# STEPS 16-18: SENSITIVITY, PHARMACOGENOMICS, SUBTYPING (depend on DE + ancestry)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 16] Submitting sensitivity analysis..."

# Build dependency: needs DE + ancestry
SENS_DEPS_LIST=()
if [[ "${JOB_DE}" != "SKIP" ]]; then SENS_DEPS_LIST+=("${JOB_DE}"); fi
if [[ "${JOB_ANCESTRY}" != "SKIP" ]]; then SENS_DEPS_LIST+=("${JOB_ANCESTRY}"); fi

SENS_DEPS=""
if [[ ${#SENS_DEPS_LIST[@]} -gt 0 ]]; then
    SENS_DEPS="afterok:$(IFS=:; echo "${SENS_DEPS_LIST[*]}")"
fi

JOB_SENSITIVITY=$(submit_job "${SCRIPT_DIR}/16_sensitivity/run_sensitivity.sbatch" "${SENS_DEPS}")
echo "  Sensitivity Job ID: ${JOB_SENSITIVITY}"

echo ""
echo "[STEP 17] Submitting pharmacogenomics..."
JOB_PHARMA=$(submit_job "${SCRIPT_DIR}/17_pharmacogenomics/run_pharma.sbatch" "${SENS_DEPS}")
echo "  Pharmacogenomics Job ID: ${JOB_PHARMA}"

echo ""
echo "[STEP 18] Submitting molecular subtyping..."
# Molecular subtyping script may not exist yet
SUBTYPING_SCRIPT="${SCRIPT_DIR}/18_molecular_subtyping/run_subtyping.sbatch"
if [[ -f "${SUBTYPING_SCRIPT}" ]]; then
    JOB_SUBTYPING=$(submit_job "${SUBTYPING_SCRIPT}" "${SENS_DEPS}")
else
    echo "  WARNING: ${SUBTYPING_SCRIPT} not found, skipping."
    JOB_SUBTYPING="SKIP"
fi
echo "  Subtyping Job ID: ${JOB_SUBTYPING}"

# ============================================================================
# STEP 19: VALIDATION (depends on ALL previous steps)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 19] Submitting validation..."

# Collect all non-SKIP job IDs for the final dependency
ALL_JOBS=()
for jid in "${JOB_QC}" "${JOB_COUNTING}" "${JOB_ANCESTRY}" "${JOB_SOMALIER}" \
           "${JOB_VARIANTS}" "${JOB_HLA}" "${JOB_DE}" "${JOB_RMATS}" \
           "${JOB_LEAFCUTTER}" "${JOB_SPLADDER}" "${JOB_BISBEE}" "${JOB_FUSIONS}" \
           "${JOB_CNV}" "${JOB_WGCNA}" "${JOB_PATHWAY}" "${JOB_IMMUNE}" \
           "${JOB_NEOANTIGEN}" "${JOB_TCR}" "${JOB_SENSITIVITY}" "${JOB_PHARMA}" \
           "${JOB_SUBTYPING}"; do
    if [[ "${jid}" != "SKIP" && -n "${jid}" ]]; then
        ALL_JOBS+=("${jid}")
    fi
done

VALIDATION_DEPS=""
if [[ ${#ALL_JOBS[@]} -gt 0 ]]; then
    VALIDATION_DEPS="afterok:$(IFS=:; echo "${ALL_JOBS[*]}")"
fi

JOB_VALIDATION=$(submit_job "${SCRIPT_DIR}/19_validation/run_validation.sbatch" "${VALIDATION_DEPS}")
echo "  Validation Job ID: ${JOB_VALIDATION}"

# ============================================================================
# STEP 20: VISUALIZATION (depends on validation)
# ============================================================================

echo "----------------------------------------------------------------------"
echo "[STEP 20] Submitting visualization..."

VIZ_DEPS=""
if [[ "${JOB_VALIDATION}" != "SKIP" ]]; then
    VIZ_DEPS="afterok:${JOB_VALIDATION}"
fi

JOB_VIZ=$(submit_job "${SCRIPT_DIR}/20_visualization/run_visualization.sbatch" "${VIZ_DEPS}")
echo "  Visualization Job ID: ${JOB_VIZ}"

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "============================================================================"
echo "  ALL JOBS SUBMITTED SUCCESSFULLY"
echo "============================================================================"
echo ""
printf "  %-35s %s\n" "Step" "Job ID"
printf "  %-35s %s\n" "-----------------------------------" "----------"
printf "  %-35s %s\n" "01 Genome Detection" "${JOB_GENOME}"
printf "  %-35s %s\n" "02 QC (array)" "${JOB_QC}"
printf "  %-35s %s\n" "03 Counting (array)" "${JOB_COUNT_ARRAY}"
printf "  %-35s %s\n" "03 Counting (merge)" "${JOB_COUNT_MERGE}"
printf "  %-35s %s\n" "04 Ancestry - GRAF extract (array)" "${JOB_GRAF_EXTRACT}"
printf "  %-35s %s\n" "04 Ancestry - GRAF merge" "${JOB_ANCESTRY}"
printf "  %-35s %s\n" "04 Ancestry - Somalier" "${JOB_SOMALIER}"
printf "  %-35s %s\n" "05 Variant Calling (array)" "${JOB_VARIANTS}"
printf "  %-35s %s\n" "06 HLA Typing (array)" "${JOB_HLA}"
printf "  %-35s %s\n" "07 Differential Expression" "${JOB_DE}"
printf "  %-35s %s\n" "08 Splicing - rMATS" "${JOB_RMATS}"
printf "  %-35s %s\n" "08 Splicing - Leafcutter" "${JOB_LEAFCUTTER}"
printf "  %-35s %s\n" "08 Splicing - SpLaDDer" "${JOB_SPLADDER}"
printf "  %-35s %s\n" "08 Splicing - BISBEE" "${JOB_BISBEE}"
printf "  %-35s %s\n" "09 Fusions (array)" "${JOB_FUSIONS}"
printf "  %-35s %s\n" "10 CNV Inference" "${JOB_CNV}"
printf "  %-35s %s\n" "11 WGCNA" "${JOB_WGCNA}"
printf "  %-35s %s\n" "12 Pathway Enrichment" "${JOB_PATHWAY}"
printf "  %-35s %s\n" "13 Immune Analysis" "${JOB_IMMUNE}"
printf "  %-35s %s\n" "14 Neoantigen Prediction" "${JOB_NEOANTIGEN}"
printf "  %-35s %s\n" "15 TCR Repertoire (array)" "${JOB_TCR}"
printf "  %-35s %s\n" "16 Sensitivity" "${JOB_SENSITIVITY}"
printf "  %-35s %s\n" "17 Pharmacogenomics" "${JOB_PHARMA}"
printf "  %-35s %s\n" "18 Molecular Subtyping" "${JOB_SUBTYPING}"
printf "  %-35s %s\n" "19 Validation" "${JOB_VALIDATION}"
printf "  %-35s %s\n" "20 Visualization" "${JOB_VIZ}"
echo ""
echo "============================================================================"
echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Cancel all:    scancel ${ALL_JOBS[*]} ${JOB_VALIDATION} ${JOB_VIZ} 2>/dev/null"
echo ""
echo "============================================================================"
