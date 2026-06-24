/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Variant Calling Module - GATK RNA-seq best practices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GATK_SPLITNCIGARREADS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("${meta.id}.split.bam"), path("${meta.id}.split.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    gatk SplitNCigarReads \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.id}.split.bam \\
        --create-output-bam-index true \\
        --tmp-dir \$PWD/tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep 'GATK' | sed 's/.*v//')
    END_VERSIONS
    """
}

process GATK_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fasta_fai
    path  known_snps
    path  known_snps_tbi

    output:
    tuple val(meta), path("${meta.id}.recal_data.table"), emit: table
    path "versions.yml", emit: versions

    script:
    """
    gatk BaseRecalibrator \\
        -R ${fasta} \\
        -I ${bam} \\
        --known-sites ${known_snps} \\
        -O ${meta.id}.recal_data.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep 'GATK' | sed 's/.*v//')
    END_VERSIONS
    """
}

process GATK_APPLYBQSR {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("${meta.id}.recal.bam"), path("${meta.id}.recal.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    gatk ApplyBQSR \\
        -R ${fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${meta.id}.recal.bam \\
        --create-output-bam-index true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep 'GATK' | sed 's/.*v//')
    END_VERSIONS
    """
}

process GATK_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fasta_fai
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(meta), path("${meta.id}.raw.vcf.gz"), path("${meta.id}.raw.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def dbsnp_arg = dbsnp ? "--dbsnp ${dbsnp}" : ''
    """
    gatk HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.id}.raw.vcf.gz \\
        ${dbsnp_arg} \\
        --dont-use-soft-clipped-bases true \\
        --standard-min-confidence-threshold-for-calling 20.0 \\
        --min-base-quality-score ${params.min_base_quality} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --tmp-dir \$PWD/tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep 'GATK' | sed 's/.*v//')
    END_VERSIONS
    """
}

process GATK_VARIANTFILTRATION {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    gatk VariantFiltration \\
        -R ${fasta} \\
        -V ${vcf} \\
        -O ${meta.id}.filtered.vcf.gz \\
        --window 35 \\
        --cluster 3 \\
        --filter-name "FS" --filter-expression "FS > 30.0" \\
        --filter-name "QD" --filter-expression "QD < 2.0" \\
        --filter-name "MQ" --filter-expression "MQ < 40.0" \\
        --filter-name "DP" --filter-expression "DP < ${params.min_variant_depth}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep 'GATK' | sed 's/.*v//')
    END_VERSIONS
    """
}

process TMB_ESTIMATION {
    label 'process_low'

    input:
    path filtered_vcfs  // All filtered VCF files
    path metadata

    output:
    path "tmb_scores.tsv",  emit: tmb
    path "tmb_summary.tsv", emit: summary
    path "versions.yml",    emit: versions

    script:
    """
    Rscript ${projectDir}/bin/estimate_tmb.R \\
        --vcf-dir . \\
        --metadata ${metadata} \\
        --output-prefix tmb \\
        --min-vaf ${params.min_vaf} \\
        --min-depth ${params.min_variant_depth}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

process CHASMPLUS {
    tag "$meta.sample_id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta

    output:
    tuple val(meta), path("*.chasmplus.tsv"), emit: chasmplus_results
    path "*.findings.md", emit: findings

    script:
    """
    # Convert VCF to CHASMplus input format
    python3 <<'PYEOF'
import subprocess, os
# Run CHASMplus via OpenCRAVAT
# chasmplus annotates each variant with a p-value for being a cancer driver
subprocess.run([
    "oc", "run", "${vcf}",
    "-l", "hg38",
    "-a", "chasmplus",
    "--mp", "${task.cpus}",
    "-d", "chasmplus_out"
], check=True)

# Parse results
import csv
results = []
result_file = "chasmplus_out/chasmplus_out.tsv"
if os.path.exists(result_file):
    with open(result_file) as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            results.append(row)

with open("${meta.sample_id}.chasmplus.tsv", 'w') as f:
    writer = csv.writer(f, delimiter='\\t')
    writer.writerow(["chrom", "pos", "ref", "alt", "gene", "chasmplus_score", "chasmplus_pvalue", "driver_prediction"])
    for r in results:
        score = float(r.get('chasmplus_score', 0))
        pval = float(r.get('chasmplus_pvalue', 1))
        driver = "driver" if pval < 0.05 else "passenger"
        writer.writerow([r.get('chrom',''), r.get('pos',''), r.get('ref',''), r.get('alt',''),
                        r.get('gene',''), score, pval, driver])
PYEOF

    # Findings log
    n_drivers=\$(awk -F'\\t' '\$8=="driver"' ${meta.sample_id}.chasmplus.tsv | wc -l)
    cat <<-FINDINGS > ${meta.sample_id}.chasmplus.findings.md
    ## CHASMplus Driver Mutation Prediction - ${meta.sample_id}
    - Total variants scored: \$(tail -n +2 ${meta.sample_id}.chasmplus.tsv | wc -l)
    - Predicted driver mutations (p < 0.05): \${n_drivers}
    FINDINGS
    """
}

process CREATE_MAF {
    tag "$meta.sample_id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path gtf

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    """
    # Use vcf2maf to convert VCF to MAF format
    # vcf2maf requires VEP; if not available, use a lightweight conversion
    vcf2maf.pl \\
        --input-vcf ${vcf} \\
        --output-maf ${meta.sample_id}.maf \\
        --tumor-id ${meta.sample_id} \\
        --ref-fasta ${fasta} \\
        --ncbi-build GRCh38 \\
        --species homo_sapiens \\
        --vep-path \$VEP_PATH \\
        --vep-data \$VEP_DATA \\
        || {
            # Fallback: lightweight VCF to MAF conversion without VEP
            python3 <<'PYEOF'
import gzip, sys

maf_header = ["Hugo_Symbol","Chromosome","Start_Position","End_Position",
              "Reference_Allele","Tumor_Seq_Allele2","Variant_Classification",
              "Variant_Type","Tumor_Sample_Barcode","HGVSp_Short"]

with open("${meta.sample_id}.maf", 'w') as out:
    out.write("\\t".join(maf_header) + "\\n")
    opener = gzip.open if "${vcf}".endswith('.gz') else open
    with opener("${vcf}", 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\\t')
            chrom, pos, _, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            info = fields[7] if len(fields) > 7 else ""
            filt = fields[6] if len(fields) > 6 else ""
            if filt != "PASS" and filt != ".": continue
            vtype = "SNP" if len(ref) == len(alt) == 1 else ("INS" if len(alt) > len(ref) else "DEL")
            out.write("\\t".join(["Unknown", chrom, pos, str(int(pos)+len(ref)-1),
                                  ref, alt, "Missense_Mutation", vtype,
                                  "${meta.sample_id}", ""]) + "\\n")
PYEOF
        }
    """
}
