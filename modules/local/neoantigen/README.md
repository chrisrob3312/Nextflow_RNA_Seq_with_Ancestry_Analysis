# Neoantigen Module

## Purpose

RNA-seq-specific neoantigen prediction from three mutation sources: SNVs/indels (pVACseq), gene fusions (NeoFuse), and alternative splicing events (SNAF). No WGS data is required. Predicts MHC class I and class II binding using patient-specific HLA types.

## Processes

| Process | Description |
|---------|-------------|
| `PVACSEQ` | Predict neoantigens from somatic SNVs/indels using pVACseq; evaluates MHC class I (8-11mer) and class II (15mer) binding via MHCflurry, NetMHCpan, NetMHCIIpan |
| `NEOFUSE` | Predict fusion-derived neoantigens from Arriba output using patient HLA alleles |
| `SNAF_SPLICING_NEOANTIGENS` | Identify neoantigens from alternative splicing junctions; predicts both T-cell and B-cell epitopes from novel splice-derived peptides |
| `MERGE_NEOANTIGENS` | Merge and deduplicate candidates from all three sources; compute per-sample neoantigen burden |

## Software

| Tool | Version | Documentation |
|------|---------|---------------|
| pVACtools (pVACseq) | 4.2.0 | [pvactools.readthedocs.io](https://pvactools.readthedocs.io/) |
| NeoFuse | 1.0 | [github.com/icbi-lab/NeoFuse](https://github.com/icbi-lab/NeoFuse) |
| SNAF | -- | [github.com/frankligy/SNAF](https://github.com/frankligy/SNAF) |

## Inputs

- `tuple val(meta), path(vcf), path(tbi)` -- Filtered variant calls (pVACseq)
- `tuple val(meta), path(fusions)` -- Arriba fusion calls (NeoFuse)
- `tuple val(meta), path(bam), path(bai)` -- Aligned BAM (SNAF)
- `path hla_json` -- HLA genotype from arcasHLA
- `path gtf` / `path snaf_db` -- Annotation and SNAF reference database

## Outputs

- `*_pvacseq/` -- pVACseq results with filtered MHC class I/II candidates
- `*_neofuse/` -- NeoFuse fusion neoantigen predictions
- `*_snaf/` -- SNAF T-cell and B-cell epitope candidates
- `merged_neoantigens.tsv` -- Combined neoantigen list from all sources
- `neoantigen_burden.tsv` -- Per-sample neoantigen burden counts

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `neoantigen_tool` | `all` | Tools: `pvacseq`, `neofuse`, `snaf`, `both`, or `all` |
| `pvac_algorithms` | `MHCflurry,NetMHCpan,NetMHCIIpan` | Binding prediction algorithms |
| `binding_threshold` | `500` | Binding affinity threshold (nM) |
| `percentile_threshold` | `2.0` | Percentile rank threshold for binding |
| `snaf_db` | `null` | Path to SNAF reference database |
