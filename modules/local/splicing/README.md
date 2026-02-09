# Splicing Module

## Purpose

Differential and aberrant splicing analysis using four complementary tools: rMATS (event-level), Leafcutter (intron cluster-level), SplAdder (graph-based event detection and testing), and Bisbee (beta-binomial modeling with protein effect prediction and outlier detection).

## Processes

| Process | Description |
|---------|-------------|
| `RMATS` | Detect differential splicing events (SE, A3SS, A5SS, MXE, RI) between two groups using a likelihood-ratio test |
| `LEAFCUTTER_JUNCTIONS` | Extract splice junctions from BAM using regtools |
| `LEAFCUTTER_CLUSTER` | Cluster introns across samples into alternatively spliced groups |
| `LEAFCUTTER_DIFF_SPLICING` | Dirichlet-multinomial GLM test for differential intron usage between groups |
| `SPLADDER_BUILD` | Construct splicing graphs from BAMs, detect and quantify AS events |
| `SPLADDER_TEST` | Negative binomial GLM differential testing of SplAdder events |
| `BISBEE_PREP` | Prepare SplAdder HDF5 junction counts for Bisbee input |
| `BISBEE_DIFF` | Beta-binomial differential splicing accounting for overdispersion |
| `BISBEE_PROT` | Predict protein effects of differential splice events (NMD, frameshifts, domain disruptions) |
| `BISBEE_OUTLIER` | Detect per-sample splicing outliers via beta-binomial model |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| rMATS | 4.3.0 | [github.com/Xinglab/rmats-turbo](https://github.com/Xinglab/rmats-turbo) |
| Leafcutter | 0.2.9 | [github.com/davidaknowles/leafcutter](https://github.com/davidaknowles/leafcutter) |
| regtools | -- | [github.com/griffithlab/regtools](https://github.com/griffithlab/regtools) |
| SplAdder | 3.0.4 | [github.com/ratschlab/spladder](https://github.com/ratschlab/spladder) |
| Bisbee | -- | [github.com/tgen/bisbee](https://github.com/tgen/bisbee) |

## Inputs

- `tuple val(meta), path(bam), path(bai)` -- Sorted BAMs (Leafcutter junctions, SplAdder)
- `path bam_list_1` / `path bam_list_2` -- Group BAM lists (rMATS)
- `path gtf` -- Gene annotation; `path fasta` -- Reference genome (Bisbee protein effects)
- `path condition_file` -- Two-column sample-to-group mapping (SplAdder test, Bisbee diff)

## Outputs

- `rmats_*/` -- Per-event-type result tables (SE, A3SS, A5SS, MXE, RI)
- `leafcutter_ds_*/` -- Cluster significance and effect sizes
- `spladder_out/` / `spladder_test_*/` -- Splicing graphs and differential test results
- `bisbee_diff_*/` -- Beta-binomial differential splicing results
- `bisbee_protein_effects/` -- NMD targets, frameshifts, domain disruptions
- `bisbee_outliers/` -- Per-sample splicing outlier calls

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `splicing_tool` | `all` | Tools to run: rmats, leafcutter, spladder, bisbee, both, or all |
| `rmats_read_length` | `150` | Read length for rMATS |
| `rmats_novel_ss` | `true` | Allow novel splice sites in rMATS |
| `leafcutter_min_coverage` | `20` | Minimum junction read support for Leafcutter clustering |
| `spladder_confidence` | `3` | SplAdder confidence level (1-3; higher = stricter) |
| `bisbee_outlier_fdr` | `0.05` | FDR threshold for Bisbee outlier detection |
