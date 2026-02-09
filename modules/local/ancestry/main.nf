/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ancestry Module - Genetic ancestry inference from bulk RNA-seq
    Tools: GRAF-anc (282K ancestry-informative SNPs), Somalier (QC/relatedness)

    Workflow:
      1. Extract genotypes at GRAF-anc 282K SNP positions from RNA-seq BAMs
      2. Merge per-sample VCFs into multi-sample PLINK format
      3. Run GRAF-anc for continental (8 groups) + subcontinental (38 groups) ancestry
      4. Also run Somalier for sample QC / relatedness checking
      5. Generate ancestry proportions (continuous) and categories (discrete)

    Note: RNA-seq typically covers 20-50K of the 282K SNPs (those in expressed
    regions). This exceeds the 10K minimum for continental inference but may be
    borderline for subcontinental resolution.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EXTRACT_GRAF_SNPS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path  graf_snp_positions  // GRAF-anc AncSnpPopAFs.txt or BED of 282K SNPs
    path  fasta

    output:
    tuple val(meta), path("${meta.id}.graf_snps.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.graf_snps.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("${meta.id}.graf_allele_counts.tsv"), emit: allele_counts
    path  "${meta.id}.snp_coverage_stats.txt",                emit: coverage_stats
    path  "versions.yml",                                      emit: versions

    script:
    """
    # Extract genotypes at GRAF-anc 282K ancestry-informative SNP positions
    # from RNA-seq BAM. Only SNPs with sufficient coverage will be called.
    # Typical yield from RNA-seq: 20,000-50,000 SNPs (in expressed regions)

    # Call variants at GRAF SNP positions with strict quality filters
    bcftools mpileup \\
        -T ${graf_snp_positions} \\
        -f ${fasta} \\
        -q 20 \\
        -Q 20 \\
        -a FORMAT/DP,FORMAT/AD \\
        --max-depth 10000 \\
        ${bam} | \\
    bcftools call \\
        -m \\
        -Oz \\
        -o ${meta.id}.graf_snps.vcf.gz

    tabix -p vcf ${meta.id}.graf_snps.vcf.gz

    # Extract allele counts for each SNP (for custom ancestry inference)
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%DP\\t%AD]\\n' \\
        ${meta.id}.graf_snps.vcf.gz > ${meta.id}.graf_allele_counts.tsv

    # Report SNP coverage statistics
    TOTAL_SNPS=\$(zcat ${graf_snp_positions} 2>/dev/null | wc -l || wc -l < ${graf_snp_positions})
    CALLED_SNPS=\$(bcftools view -H ${meta.id}.graf_snps.vcf.gz | wc -l)
    HQ_SNPS=\$(bcftools view -H -i 'FORMAT/DP>=10' ${meta.id}.graf_snps.vcf.gz | wc -l)
    echo -e "sample_id\\ttotal_graf_snps\\tcalled_snps\\thq_snps_dp10\\tpct_called" > ${meta.id}.snp_coverage_stats.txt
    PCT=\$(echo "scale=1; \${CALLED_SNPS} * 100 / \${TOTAL_SNPS}" | bc 2>/dev/null || echo "NA")
    echo -e "${meta.id}\\t\${TOTAL_SNPS}\\t\${CALLED_SNPS}\\t\${HQ_SNPS}\\t\${PCT}" >> ${meta.id}.snp_coverage_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process MERGE_GRAF_VCFS {
    label 'process_medium'

    input:
    path vcf_files   // All per-sample GRAF SNP VCFs
    path tbi_files   // All indices

    output:
    path "all_samples_graf_snps.vcf.gz",     emit: merged_vcf
    path "all_samples_graf_snps.vcf.gz.tbi", emit: merged_tbi
    path "versions.yml",                      emit: versions

    script:
    """
    # Merge all per-sample VCFs into a multi-sample VCF
    bcftools merge \\
        --force-samples \\
        -Oz \\
        -o all_samples_graf_snps.vcf.gz \\
        *.graf_snps.vcf.gz

    tabix -p vcf all_samples_graf_snps.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process GRAFANC_RUN {
    label 'process_medium'

    input:
    path merged_vcf
    path merged_tbi
    path grafanc_data  // GRAF-anc data directory (AncSnpPopAFs.txt, etc.)

    output:
    path "grafanc_results.txt",             emit: results
    path "grafanc_ancestry_summary.tsv",    emit: summary
    path "versions.yml",                    emit: versions

    script:
    """
    # Run GRAF-anc on the merged multi-sample VCF
    # GRAF-anc outputs: sample ID, SNP count, GD1-GD3 (continental distances),
    # subcontinental scores (EA1-4, AF1-3, EU1-3, SA1-2, IC1-3),
    # ancestry proportions (Pe, Pf, Pa), and AncGroupID

    grafanc \\
        -vcf ${merged_vcf} \\
        -data ${grafanc_data} \\
        -out grafanc_results.txt \\
        -thread ${task.cpus}

    # Parse GRAF-anc output into standardized format
    python3 ${projectDir}/bin/parse_grafanc_output.py \\
        --input grafanc_results.txt \\
        --output grafanc_ancestry_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grafanc: \$(grafanc 2>&1 | grep -i version | head -1 || echo 'unknown')
    END_VERSIONS
    """
}

process ANCESTRY_INFERENCE {
    label 'process_medium'

    input:
    path grafanc_results    // GRAF-anc output (if available)
    path allele_counts      // All sample allele count files (fallback)
    path reference_panel    // Reference panel with known ancestry
    path reference_labels   // Population labels

    output:
    path "ancestry_proportions.tsv",  emit: proportions
    path "ancestry_categories.tsv",   emit: categories
    path "ancestry_pca.tsv",          emit: pca
    path "ancestry_plots/",           emit: plots
    path "ancestry_findings.md",      emit: findings
    path "versions.yml",              emit: versions

    script:
    """
    python3 ${projectDir}/bin/ancestry_inference.py \\
        --grafanc-results ${grafanc_results} \\
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

process SOMALIER_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  sites_vcf   // Somalier sites VCF (~17K coding-region SNPs)

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

process SOMALIER_RELATE {
    label 'process_low'

    input:
    path extracted_files  // All .somalier files

    output:
    path "somalier.samples.tsv", emit: samples
    path "somalier.pairs.tsv",   emit: pairs
    path "somalier.html",        emit: html
    path "versions.yml",         emit: versions

    script:
    """
    somalier relate \\
        -o somalier \\
        ${extracted_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier 2>&1 | grep -o 'version [0-9.]*' | sed 's/version //')
    END_VERSIONS
    """
}
