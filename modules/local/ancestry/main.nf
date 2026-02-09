/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ancestry Module - Genetic ancestry inference from bulk RNA-seq
    Tools: Somalier (extract/relate/ancestry), GRAF-pop SNP extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SOMALIER_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  sites_vcf   // Somalier sites VCF (ancestry-informative SNPs)

    output:
    tuple val(meta), path("*.somalier"), emit: extracted
    path "versions.yml",                 emit: versions

    script:
    """
    somalier extract \\
        --sites ${sites_vcf} \\
        -f ${fasta} \\
        -d . \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier 2>&1 | grep -o 'version [0-9.]*' | sed 's/version //')
    END_VERSIONS
    """
}

process SOMALIER_ANCESTRY {
    label 'process_low'

    input:
    path extracted_files  // All .somalier files
    path labels           // Ancestry labels for reference panel
    path reference_somalier  // Reference panel .somalier files directory

    output:
    path "somalier-ancestry.somalier-ancestry.tsv", emit: ancestry_tsv
    path "somalier-ancestry*.html",                 emit: html
    path "versions.yml",                            emit: versions

    script:
    """
    somalier ancestry \\
        --labels ${labels} \\
        ${reference_somalier}/*.somalier ++ \\
        ${extracted_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier 2>&1 | grep -o 'version [0-9.]*' | sed 's/version //')
    END_VERSIONS
    """
}

process EXTRACT_GRAF_SNPS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path  graf_snp_bed
    path  fasta

    output:
    tuple val(meta), path("${meta.id}.graf_snps.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.graf_snps.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("${meta.id}.graf_allele_counts.tsv"), emit: allele_counts
    path "versions.yml",                                       emit: versions

    script:
    """
    # Extract reads at GRAF ancestry-informative SNP positions
    samtools mpileup \\
        -l ${graf_snp_bed} \\
        -f ${fasta} \\
        -q 20 \\
        -Q 20 \\
        --output-tags DP,AD \\
        ${bam} | \\
    bcftools call -m -Oz -o ${meta.id}.graf_snps.vcf.gz

    tabix -p vcf ${meta.id}.graf_snps.vcf.gz

    # Extract allele counts at each GRAF SNP for ancestry inference
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%DP\\t%AD]\\n' \\
        ${meta.id}.graf_snps.vcf.gz > ${meta.id}.graf_allele_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process ANCESTRY_INFERENCE {
    label 'process_medium'

    input:
    path allele_counts    // All sample allele count files
    path reference_panel  // Reference panel with known ancestry
    path reference_labels // Population labels

    output:
    path "ancestry_proportions.tsv",  emit: proportions
    path "ancestry_categories.tsv",   emit: categories
    path "ancestry_pca.tsv",          emit: pca
    path "ancestry_plots/",           emit: plots
    path "versions.yml",              emit: versions

    script:
    """
    python3 ${projectDir}/bin/ancestry_inference.py \\
        --allele-counts ${allele_counts} \\
        --reference-panel ${reference_panel} \\
        --reference-labels ${reference_labels} \\
        --categories ${params.graf_ancestry_categories} \\
        --output-prefix ancestry \\
        --plot-dir ancestry_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}
