/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Counting Module - featureCounts (Subread) and normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SUBREAD_FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path  gtf

    output:
    tuple val(meta), path("${meta.id}.featureCounts.txt"),         emit: counts
    tuple val(meta), path("${meta.id}.featureCounts.txt.summary"), emit: summary
    path "versions.yml",                                           emit: versions

    script:
    def strandedness = params.fc_strandedness
    def extra_attr   = params.fc_extra_attributes ? "-g ${params.fc_extra_attributes}" : ''
    def pe_flag      = '-p --countReadPairs'
    """
    featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o ${meta.id}.featureCounts.txt \\
        -s ${strandedness} \\
        -t ${params.fc_count_type} \\
        -g ${params.fc_group_features} \\
        ${extra_attr} \\
        ${pe_flag} \\
        -B \\
        -C \\
        --largestOverlap \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -o 'v[0-9.]*' | sed 's/v//')
    END_VERSIONS
    """
}

process MERGE_COUNTS {
    label 'process_low'

    input:
    path count_files  // All individual featureCounts output files

    output:
    path "raw_count_matrix.tsv",    emit: count_matrix
    path "gene_annotations.tsv",    emit: gene_info
    path "count_summary_stats.tsv", emit: summary_stats
    path "versions.yml",            emit: versions

    script:
    """
    Rscript ${projectDir}/bin/merge_counts.R \\
        --input-dir . \\
        --pattern "*.featureCounts.txt" \\
        --output-matrix raw_count_matrix.tsv \\
        --output-genes gene_annotations.tsv \\
        --output-stats count_summary_stats.tsv \\
        --min-counts ${params.min_gene_counts} \\
        --min-samples ${params.min_samples_expressing}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

process NORMALIZE_COUNTS {
    label 'process_medium'

    input:
    path count_matrix
    path metadata

    output:
    path "normalized_counts_${params.normalization_method}.tsv", emit: normalized
    path "normalized_counts_tpm.tsv",                            emit: tpm
    path "normalization_qc/",                                    emit: qc_plots
    path "versions.yml",                                         emit: versions

    script:
    """
    Rscript ${projectDir}/bin/normalize_counts.R \\
        --counts ${count_matrix} \\
        --metadata ${metadata} \\
        --method ${params.normalization_method} \\
        --output-prefix normalized_counts \\
        --plot-dir normalization_qc \\
        --min-counts ${params.min_gene_counts} \\
        --min-samples ${params.min_samples_expressing}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        DESeq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}

process BATCH_CORRECTION {
    label 'process_medium'

    input:
    path count_matrix
    path metadata

    output:
    path "batch_corrected_counts.tsv", emit: corrected_counts
    path "batch_correction_qc/",       emit: qc_plots
    path "versions.yml",               emit: versions

    script:
    """
    Rscript ${projectDir}/bin/batch_correction.R \\
        --counts ${count_matrix} \\
        --metadata ${metadata} \\
        --method ${params.batch_correction} \\
        --batch-variable ${params.batch_variable} \\
        --output batch_corrected_counts.tsv \\
        --plot-dir batch_correction_qc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        sva: \$(Rscript -e "cat(as.character(packageVersion('sva')))")
    END_VERSIONS
    """
}
