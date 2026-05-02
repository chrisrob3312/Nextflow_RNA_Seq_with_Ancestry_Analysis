# SBATCH Version - Cancer RNA-Seq Pipeline

Standalone SLURM batch job alternative to the Nextflow cancer RNA-seq pipeline. This version submits individual sbatch scripts with explicit dependency chains, giving you full control over each step on an HPC cluster without requiring Nextflow or container runtimes.

## Configuration

Edit `config.sh` before running the pipeline:

1. Set `PROJECT_DIR` to your project root directory
2. Update all reference file paths (`FASTA`, `GTF`, `KNOWN_SNPS`, etc.)
3. Set `SAMPLESHEET` to point to your sample CSV (format: sample_id, bam_path, bai_path, ...)
4. Configure SLURM settings (`SLURM_ACCOUNT`, `SLURM_PARTITION`)
5. Adjust analysis parameters as needed (thresholds, covariates, tools)
6. Uncomment and configure your module/conda environment loading approach

## How to Run

### Full pipeline

```bash
# 1. Edit config.sh with your paths and settings
# 2. Run the master submission script:
./submit_all.sh
```

This will validate inputs, then submit all steps with SLURM `--dependency=afterok` flags so that each step waits for its prerequisites.

### Individual steps

Submit any single step manually:

```bash
# Example: just run QC
sbatch --array=0-$((N_SAMPLES-1)) 02_qc/run_qc.sbatch

# Example: run DE after counting is done (provide the counting job ID)
sbatch --dependency=afterok:12345678 07_differential_expression/run_de.sbatch
```

## Directory Structure

```
SBATCH_Version/
├── config.sh                      # Shared configuration (edit this first)
├── submit_all.sh                  # Master submission with dependency chains
├── README.md                      # This file
├── 00_setup/                      # Directory creation and input validation
├── 01_genome_detection/           # Detect hg19 vs hg38 from BAM headers
├── 02_qc/                         # Samtools, RSeQC, Picard, MultiQC
├── 03_counting/                   # featureCounts + normalization
├── 04_ancestry/                   # GRAF-pop + Somalier ancestry inference
├── 05_variant_calling/            # GATK RNA-seq variant calling
├── 06_hla_typing/                 # HLA allele typing from RNA-seq
├── 07_differential_expression/    # DESeq2 + limma-voom
├── 08_splicing/                   # rMATS, Leafcutter, SpLaDDer, BISBEE
├── 09_fusions/                    # FusionCatcher + Arriba
├── 10_cnv_inference/              # inferCNV from expression
├── 11_wgcna/                      # Weighted gene co-expression networks
├── 12_pathway_enrichment/         # GO, KEGG, Reactome, Hallmark, etc.
├── 13_immune_analysis/            # Immune deconvolution (multi-method)
├── 14_neoantigen/                 # Neoantigen prediction
├── 15_tcr_repertoire/             # TRUST4 TCR/BCR extraction
├── 16_sensitivity/                # Drug sensitivity prediction
├── 17_pharmacogenomics/           # Pharmacogenomic variant annotation
├── 18_molecular_subtyping/        # Molecular subtype classification
├── 19_validation/                 # Cross-method concordance + QC checks
└── 20_visualization/              # Integrated figures + HTML index
```

## Re-running a Failed Step

1. Check the SLURM log for the failed job:
   ```bash
   cat logs/<step_name>/<script>_<jobid>.err
   ```

2. Fix the issue (missing input, resource limits, software error).

3. Re-submit just that step with appropriate dependencies:
   ```bash
   # If the step has no unmet dependencies (upstream completed):
   sbatch 07_differential_expression/run_de.sbatch

   # If it is an array job and only some tasks failed:
   sbatch --array=3,7,12 05_variant_calling/gatk_rnaseq.sbatch
   ```

4. Re-submit downstream steps with `--dependency=afterok:<new_job_id>`.

Note: There is no automatic resume. You must manually re-submit failed and downstream steps.

## Dependency Graph

```
setup (local)
  └── genome_detection
        ├── qc (array)
        └── counting (array → merge)
              ├── ancestry (array → merge) ──┐
              │     └─────────────────────────┼── DE
              ├── variant_calling (array) ────┼──────────────┐
              ├── hla_typing (array) ─────────┼──────────────┼── neoantigen
              ├── splicing (4 tools, parallel) │              │
              ├── fusions (array)             │              │
              ├── cnv                         │              │
              ├── immune                      │              │
              ├── tcr (array)                 │              │
              │                               │              │
              │                        DE ────┤              │
              │                          ├── wgcna           │
              │                          ├── pathway         │
              │                          │                   │
              │                   DE + ancestry ──┐          │
              │                                   ├── sensitivity
              │                                   ├── pharmacogenomics
              │                                   └── subtyping
              │
              └── [all above] ──── validation ──── visualization
```

## Differences from the Nextflow Version

| Feature | Nextflow | SBATCH Version |
|---------|----------|----------------|
| Automatic resume on failure | Yes (`-resume`) | No - re-submit manually |
| Containerization | Docker/Singularity per process | Load modules or activate conda yourself |
| Dependency management | Implicit via channels | Explicit `--dependency=afterok` flags |
| Portability | Any system with Nextflow | Requires SLURM scheduler |
| Resource tuning | Per-process in nextflow.config | Per-script in SBATCH headers |
| Parallelism | Automatic | Manual via array jobs |
| Provenance tracking | Built-in reports | Manual (check SLURM logs) |

### When to use this version

- Your HPC does not support Nextflow or containers
- You need fine-grained control over each step's resource allocation
- You want to inspect/modify intermediate outputs between steps
- You prefer explicit SLURM workflows over Nextflow DSL
