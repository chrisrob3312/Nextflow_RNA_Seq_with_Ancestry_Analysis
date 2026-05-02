#!/bin/bash
# ============================================================================
# SHARED CONFIGURATION FOR SBATCH PIPELINE
# ============================================================================
# Source this file at the top of every sbatch script:
#   source "$(dirname "$0")/../config.sh"
#
# Edit paths below to match your HPC environment.
# ============================================================================

# --- Project paths ---
export PROJECT_DIR="/path/to/project"
export BAM_DIR="${PROJECT_DIR}/bams"
export RESULTS_DIR="${PROJECT_DIR}/results"
export LOGS_DIR="${PROJECT_DIR}/logs"
export SCRIPTS_DIR="${PROJECT_DIR}/SBATCH_Version"
export TEMP_DIR="${PROJECT_DIR}/tmp"

# --- Reference files ---
export GENOME="GRCh38"
export FASTA="/path/to/refs/GRCh38.primary_assembly.genome.fa"
export FASTA_FAI="${FASTA}.fai"
export GTF="/path/to/refs/gencode.v44.primary_assembly.annotation.gtf"
export GENE_BED="/path/to/refs/gencode.v44.genes.bed12"
export KNOWN_SNPS="/path/to/refs/dbsnp_146.hg38.vcf.gz"
export KNOWN_SNPS_TBI="${KNOWN_SNPS}.tbi"
export CHAIN_FILE="/path/to/refs/hg19ToHg38.over.chain.gz"

# --- Ancestry references ---
export GRAF_SNP_BED="/path/to/refs/GrafAncSnpFile.txt"
export GRAFANC_DATA="/path/to/refs/grafanc_data"
export GRAFANC_BINARY="/path/to/software/grafpop"
export SOMALIER_SITES="/path/to/refs/somalier_sites.hg38.vcf.gz"
export SOMALIER_ANCESTRY_LABELS="/path/to/refs/somalier/ancestry-labels-1kg.tsv"

# --- Fusion references ---
export FUSIONCATCHER_DATA="/path/to/refs/fusioncatcher_data"
export ARRIBA_BLACKLIST="/path/to/refs/arriba_blacklist_hg38_GRCh38_v2.4.0.tsv.gz"
export ARRIBA_KNOWN_FUSIONS="/path/to/refs/arriba_known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"
export ARRIBA_PROTEIN_DOMAINS="/path/to/refs/arriba_protein_domains_hg38_GRCh38_v2.4.0.gff3"

# --- Neoantigen references ---
export SNAF_DB="/path/to/refs/snaf_db"

# --- Tool paths (if not using modules) ---
# export SAMTOOLS="/path/to/samtools"
# export BCFTOOLS="/path/to/bcftools"

# --- Sample sheet ---
export SAMPLESHEET="${PROJECT_DIR}/samplesheet.csv"

# --- Analysis parameters ---
export TARGET_BUILD="hg38"
export FC_STRANDEDNESS=2              # 0=unstranded, 1=stranded, 2=reverse
export MIN_GENE_COUNTS=10
export MIN_SAMPLES_EXPRESSING=3
export NORMALIZATION_METHOD="vst"
export BATCH_VARIABLE="batch"
export DE_PADJ=0.05
export DE_LFC=0.585                   # log2(1.5)
export DE_COVARIATES="batch,sex,age,blast_percentage,tumor_purity,timepoint"
export SPLICING_TOOL="all"            # rmats, leafcutter, spladder, bisbee, all
export RMATS_READ_LENGTH=150
export LEAFCUTTER_MIN_COV=20
export SPLADDER_CONFIDENCE=3
export BISBEE_OUTLIER_FDR=0.05
export FUSION_MIN_READS=3
export HLA_GENES="A,B,C,DPB1,DQB1,DQA1,DRB1"
export BINDING_THRESHOLD=500
export WGCNA_MIN_MODULE_SIZE=30
export WGCNA_MERGE_HEIGHT=0.25
export WGCNA_NETWORK_TYPE="signed"
export BOOTSTRAP_ITERATIONS=1000
export OUTLIER_SD_THRESHOLD=3
export MIN_SAMPLE_CORRELATION=0.8
export PATHWAY_DBS="GO_BP,GO_MF,KEGG,REACTOME,HALLMARK,IMMUNESIGDB"
export DRUG_RESPONSE_DB="GDSC"

# --- SLURM defaults ---
export SLURM_ACCOUNT="your_account"
export SLURM_PARTITION="normal"
export SLURM_PARTITION_BIGMEM="bigmem"
export SLURM_EMAIL="your_email@institution.edu"

# --- Conda/module environment ---
# Uncomment the approach you use:
# export CONDA_ENV_BASE="/path/to/conda/envs"
# module use /path/to/modulefiles

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

get_sample_ids() {
    tail -n +2 "${SAMPLESHEET}" | cut -d',' -f1
}

get_bam_for_sample() {
    local sample_id="$1"
    grep "^${sample_id}," "${SAMPLESHEET}" | cut -d',' -f2
}

get_bai_for_sample() {
    local sample_id="$1"
    grep "^${sample_id}," "${SAMPLESHEET}" | cut -d',' -f3
}

get_metadata_field() {
    local sample_id="$1"
    local field_num="$2"
    grep "^${sample_id}," "${SAMPLESHEET}" | cut -d',' -f"${field_num}"
}

check_file_exists() {
    local file="$1"
    local desc="$2"
    if [[ ! -f "$file" ]]; then
        echo "ERROR: ${desc} not found: ${file}" >&2
        return 1
    fi
}

log_step() {
    local step="$1"
    local msg="$2"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${step}] ${msg}" | tee -a "${LOGS_DIR}/pipeline.log"
}

mkdir -p "${RESULTS_DIR}" "${LOGS_DIR}" "${TEMP_DIR}"
