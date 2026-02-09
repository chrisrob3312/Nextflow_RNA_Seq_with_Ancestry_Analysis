/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Neoantigen Prediction Module - pVACseq and NeoFuse (RNA-seq specific)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PVACSEQ {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(hla_json)

    output:
    tuple val(meta), path("${meta.id}_pvacseq/"), emit: results
    tuple val(meta), path("${meta.id}_pvacseq/MHC_Class_I/${meta.id}.filtered.tsv"), emit: class_i, optional: true
    tuple val(meta), path("${meta.id}_pvacseq/MHC_Class_II/${meta.id}.filtered.tsv"), emit: class_ii, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Parse HLA types from arcasHLA JSON
    HLA_ALLELES=\$(python3 ${projectDir}/bin/parse_hla_for_pvac.py ${hla_json})

    mkdir -p ${meta.id}_pvacseq

    pvacseq run \\
        ${vcf} \\
        ${meta.id} \\
        \${HLA_ALLELES} \\
        ${params.pvac_algorithms} \\
        ${meta.id}_pvacseq \\
        -e1 8,9,10,11 \\
        -e2 15 \\
        --iedb-install-directory /opt/iedb \\
        -b ${params.binding_threshold} \\
        --percentile-threshold ${params.percentile_threshold} \\
        -t ${task.cpus} \\
        --normal-sample-name normal \\
        --keep-tmp-files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvactools: \$(pvacseq --version 2>&1 | sed 's/.*pvacseq //')
    END_VERSIONS
    """
}

process NEOFUSE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(fusions)   // Arriba fusion output
    tuple val(meta2), path(hla_json)

    output:
    tuple val(meta), path("${meta.id}_neofuse/"), emit: results
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p ${meta.id}_neofuse

    # Parse HLA types
    HLA_ALLELES=\$(python3 ${projectDir}/bin/parse_hla_for_pvac.py ${hla_json})

    NeoFuse \\
        --fusions ${fusions} \\
        --hla-alleles \${HLA_ALLELES} \\
        --output ${meta.id}_neofuse \\
        --sample-id ${meta.id} \\
        --binding-threshold ${params.binding_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neofuse: \$(NeoFuse --version 2>&1 | head -1 || echo 'unknown')
    END_VERSIONS
    """
}

process MERGE_NEOANTIGENS {
    label 'process_low'

    input:
    path pvacseq_results
    path neofuse_results

    output:
    path "merged_neoantigens.tsv",      emit: merged
    path "neoantigen_summary.tsv",      emit: summary
    path "neoantigen_burden.tsv",       emit: burden
    path "versions.yml",                emit: versions

    script:
    """
    python3 ${projectDir}/bin/neoantigen_filter.py \\
        --pvacseq-dir . \\
        --neofuse-dir . \\
        --binding-threshold ${params.binding_threshold} \\
        --output-merged merged_neoantigens.tsv \\
        --output-summary neoantigen_summary.tsv \\
        --output-burden neoantigen_burden.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
