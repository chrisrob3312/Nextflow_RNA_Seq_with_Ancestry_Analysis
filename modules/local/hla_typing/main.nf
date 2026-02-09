/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HLA Typing Module - arcasHLA and OptiType from RNA-seq BAMs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ARCASHLA_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.extracted.1.fq.gz"), path("${meta.id}.extracted.2.fq.gz"), emit: extracted
    path "versions.yml", emit: versions

    script:
    """
    arcasHLA extract \\
        ${bam} \\
        -o . \\
        -t ${task.cpus} \\
        --paired \\
        -v

    # Rename outputs
    mv *.extracted.1.fq.gz ${meta.id}.extracted.1.fq.gz
    mv *.extracted.2.fq.gz ${meta.id}.extracted.2.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcasHLA: \$(arcasHLA --version 2>&1 | sed 's/.*arcasHLA //')
    END_VERSIONS
    """
}

process ARCASHLA_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fq1), path(fq2)

    output:
    tuple val(meta), path("${meta.id}.genotype.json"), emit: genotype
    path "versions.yml", emit: versions

    script:
    """
    arcasHLA genotype \\
        ${fq1} ${fq2} \\
        -g ${params.arcashla_genes} \\
        -o . \\
        -t ${task.cpus} \\
        -v

    mv *.genotype.json ${meta.id}.genotype.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcasHLA: \$(arcasHLA --version 2>&1 | sed 's/.*arcasHLA //')
    END_VERSIONS
    """
}

process ARCASHLA_MERGE {
    label 'process_low'

    input:
    path genotype_jsons

    output:
    path "hla_genotypes.tsv", emit: merged
    path "versions.yml",      emit: versions

    script:
    """
    arcasHLA merge \\
        --run . \\
        -o hla_genotypes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcasHLA: \$(arcasHLA --version 2>&1 | sed 's/.*arcasHLA //')
    END_VERSIONS
    """
}

process OPTITYPE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_optitype/"), emit: results
    tuple val(meta), path("${meta.id}_optitype/*_result.tsv"), emit: typing
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p ${meta.id}_optitype

    # Extract HLA reads
    samtools view -@ ${task.cpus} -h ${bam} chr6:29941260-29945884 chr6:31353872-31357187 chr6:31268749-31272092 | \\
        samtools sort -n -@ ${task.cpus} | \\
        samtools fastq -@ ${task.cpus} \\
            -1 hla_R1.fq.gz \\
            -2 hla_R2.fq.gz \\
            -0 /dev/null -s /dev/null -

    OptiTypePipeline.py \\
        -i hla_R1.fq.gz hla_R2.fq.gz \\
        --rna \\
        -o ${meta.id}_optitype \\
        -v

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        optitype: \$(OptiTypePipeline.py --version 2>&1 | sed 's/.*OptiType //')
    END_VERSIONS
    """
}
