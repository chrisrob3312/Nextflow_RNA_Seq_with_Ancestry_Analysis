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
