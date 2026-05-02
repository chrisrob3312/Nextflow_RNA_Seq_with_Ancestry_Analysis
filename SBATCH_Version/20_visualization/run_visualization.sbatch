#!/bin/bash
# ============================================================================
# INTEGRATED VISUALIZATION
# ============================================================================
# Generates summary figures and an HTML index from all pipeline outputs.
#
# Submit: sbatch run_visualization.sbatch
# ============================================================================

#SBATCH --job-name=visualization
#SBATCH --output=../logs/visualization/visualization_%j.out
#SBATCH --error=../logs/visualization/visualization_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
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
# PARAMETERS
# ============================================================================

RESULTS_INPUT="${RESULTS_INPUT:-${RESULTS_DIR}}"
OUTPUT_DIR="${VIZ_OUTPUT_DIR:-${RESULTS_DIR}/20_visualization}"

# ============================================================================
# VALIDATE INPUTS
# ============================================================================

log_step "VIZ" "Starting integrated visualization"

if [[ ! -d "${RESULTS_INPUT}" ]]; then
    echo "ERROR: Results directory not found: ${RESULTS_INPUT}"
    exit 1
fi

# ============================================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================================

mkdir -p "${OUTPUT_DIR}"/{figures,html}
mkdir -p "$(dirname "$0")/logs"

# ============================================================================
# LOG PARAMETERS
# ============================================================================

log_step "VIZ" "Parameters:"
log_step "VIZ" "  Results dir: ${RESULTS_INPUT}"
log_step "VIZ" "  Output dir: ${OUTPUT_DIR}"
log_step "VIZ" "  CPUs: ${SLURM_CPUS_PER_TASK:-4}"

# ============================================================================
# LOAD R ENVIRONMENT
# ============================================================================

# Uncomment whichever module/conda system your cluster uses:
# module load R/4.3.0
# conda activate rnaseq_viz

# ============================================================================
# RUN VISUALIZATION SCRIPT
# ============================================================================

Rscript "${SCRIPT_DIR}/create_visualizations.R" \
    --results_dir "${RESULTS_INPUT}" \
    --output_dir "${OUTPUT_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-4}"

# ============================================================================
# CHECK COMPLETION
# ============================================================================

if [[ $? -eq 0 ]]; then
    log_step "VIZ" "Visualization generation completed successfully"
    log_step "VIZ" "Figures: ${OUTPUT_DIR}/figures/"
    log_step "VIZ" "Index: ${OUTPUT_DIR}/html/figure_index.html"
else
    log_step "VIZ" "ERROR: Visualization generation failed"
    exit 1
fi

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "============================================"
echo "Visualization Summary"
echo "============================================"
N_FIGURES=$(find "${OUTPUT_DIR}/figures" -name "*.png" -o -name "*.pdf" 2>/dev/null | wc -l)
echo "  Figures generated: ${N_FIGURES}"
if [[ -f "${OUTPUT_DIR}/html/figure_index.html" ]]; then
    echo "  HTML index: ${OUTPUT_DIR}/html/figure_index.html"
fi
echo "============================================"

log_step "VIZ" "Job complete"
