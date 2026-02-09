/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Neoantigen Prediction Module - pVACseq, NeoFuse, and SNAF
    RNA-seq specific neoantigen identification (no WGS required)
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
    path "pvacseq_findings.md", emit: findings
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

    # Findings
    CLASS_I=\$(wc -l < ${meta.id}_pvacseq/MHC_Class_I/${meta.id}.filtered.tsv 2>/dev/null || echo 1)
    CLASS_II=\$(wc -l < ${meta.id}_pvacseq/MHC_Class_II/${meta.id}.filtered.tsv 2>/dev/null || echo 1)
    cat <<-FINDINGS > pvacseq_findings.md
    ### pVACseq Neoantigen Prediction: ${meta.id}
    - MHC Class I candidates: \$((CLASS_I - 1))
    - MHC Class II candidates: \$((CLASS_II - 1))
    - Binding threshold: ${params.binding_threshold} nM
    FINDINGS

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
    path "neofuse_findings.md", emit: findings
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

    # Findings
    NEOAG=\$(wc -l < ${meta.id}_neofuse/NeoFuse_results.tsv 2>/dev/null || echo 1)
    cat <<-FINDINGS > neofuse_findings.md
    ### NeoFuse Fusion Neoantigens: ${meta.id}
    - Fusion neoantigen candidates: \$((NEOAG - 1))
    FINDINGS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neofuse: \$(NeoFuse --version 2>&1 | head -1 || echo 'unknown')
    END_VERSIONS
    """
}

process SNAF_SPLICING_NEOANTIGENS {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path  gtf
    path  hla_json
    path  snaf_db   // SNAF reference database

    output:
    tuple val(meta), path("${meta.id}_snaf/"),                  emit: results
    tuple val(meta), path("${meta.id}_snaf/T_candidates.tsv"),  emit: t_cell, optional: true
    tuple val(meta), path("${meta.id}_snaf/B_candidates.tsv"),  emit: b_cell, optional: true
    path "snaf_findings.md",                                    emit: findings
    path "versions.yml",                                        emit: versions

    script:
    """
    mkdir -p ${meta.id}_snaf

    # SNAF (Splicing Neo Antigen Finder)
    # Identifies neoantigens derived from alternative splicing events
    # Predicts both T-cell and B-cell epitopes from splice junctions
    snaf analyze \\
        --bam ${bam} \\
        --gtf ${gtf} \\
        --hla ${hla_json} \\
        --db ${snaf_db} \\
        --outdir ${meta.id}_snaf \\
        --sample ${meta.id} \\
        --cores ${task.cpus}

    # Findings
    T_NEOAG=\$(wc -l < ${meta.id}_snaf/T_candidates.tsv 2>/dev/null || echo 1)
    B_NEOAG=\$(wc -l < ${meta.id}_snaf/B_candidates.tsv 2>/dev/null || echo 1)
    cat <<-FINDINGS > snaf_findings.md
    ### SNAF Splicing Neoantigens: ${meta.id}
    - T-cell epitope candidates: \$((T_NEOAG - 1))
    - B-cell epitope candidates: \$((B_NEOAG - 1))
    - Source: alternative splicing junctions in RNA-seq
    FINDINGS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snaf: \$(snaf --version 2>&1 | head -1 || echo 'unknown')
    END_VERSIONS
    """
}

process MERGE_NEOANTIGENS {
    label 'process_low'

    input:
    path pvacseq_results
    path neofuse_results
    path snaf_results

    output:
    path "merged_neoantigens.tsv",      emit: merged
    path "neoantigen_summary.tsv",      emit: summary
    path "neoantigen_burden.tsv",       emit: burden
    path "neoantigen_findings.md",      emit: findings
    path "versions.yml",                emit: versions

    script:
    """
    python3 ${projectDir}/bin/neoantigen_filter.py \\
        --pvacseq-dir . \\
        --neofuse-dir . \\
        --snaf-dir . \\
        --binding-threshold ${params.binding_threshold} \\
        --output-merged merged_neoantigens.tsv \\
        --output-summary neoantigen_summary.tsv \\
        --output-burden neoantigen_burden.tsv

    # Findings
    TOTAL=\$(awk 'NR>1' merged_neoantigens.tsv 2>/dev/null | wc -l || echo 0)
    cat <<-FINDINGS > neoantigen_findings.md
    ### Merged Neoantigen Summary
    - Total unique neoantigen candidates: \${TOTAL}
    - Sources: pVACseq (SNV/indel), NeoFuse (fusion), SNAF (splicing)
    - See neoantigen_burden.tsv for per-sample burden
    FINDINGS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
