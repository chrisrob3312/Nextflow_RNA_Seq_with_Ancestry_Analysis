/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Differential Splicing Module - rMATS and Leafcutter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RMATS {
    label 'process_high'

    input:
    path bam_list_1  // Text file listing BAM paths for group 1
    path bam_list_2  // Text file listing BAM paths for group 2
    path gtf
    val  contrast_name

    output:
    path "rmats_${contrast_name}/", emit: results
    path "versions.yml",            emit: versions

    script:
    """
    mkdir -p rmats_${contrast_name}

    rmats.py \\
        --b1 ${bam_list_1} \\
        --b2 ${bam_list_2} \\
        --gtf ${gtf} \\
        -t paired \\
        --readLength ${params.rmats_read_length} \\
        --nthread ${task.cpus} \\
        --od rmats_${contrast_name} \\
        --tmp rmats_tmp \\
        ${params.rmats_novel_ss ? '--novelSS' : ''} \\
        --statoff false

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(rmats.py --version 2>&1 | sed 's/.*v//')
    END_VERSIONS
    """
}

process LEAFCUTTER_JUNCTIONS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.junc"), emit: junctions
    path "versions.yml",                      emit: versions

    script:
    """
    regtools junctions extract \\
        -a 8 \\
        -m 50 \\
        -M 500000 \\
        -s 0 \\
        ${bam} > ${meta.id}.junc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(regtools --version 2>&1 | head -1)
    END_VERSIONS
    """
}

process LEAFCUTTER_CLUSTER {
    label 'process_medium'

    input:
    path junc_files

    output:
    path "leafcutter_perind_numers.counts.gz", emit: counts
    path "versions.yml",                        emit: versions

    script:
    """
    # Create junction file list
    ls *.junc > junc_files.txt

    python3 \$(which leafcutter_cluster_regtools.py) \\
        -j junc_files.txt \\
        -m ${params.leafcutter_min_coverage} \\
        -o leafcutter \\
        -l 500000

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leafcutter: \$(pip show leafcutter 2>/dev/null | grep Version | sed 's/Version: //' || echo 'unknown')
    END_VERSIONS
    """
}

process LEAFCUTTER_DIFF_SPLICING {
    label 'process_high'

    input:
    path counts
    path groups_file  // Two-column: sample_id, group
    val  contrast_name

    output:
    path "leafcutter_ds_${contrast_name}/", emit: results
    path "versions.yml",                    emit: versions

    script:
    """
    mkdir -p leafcutter_ds_${contrast_name}

    Rscript \$(which leafcutter_ds.R) \\
        --num_threads ${task.cpus} \\
        ${counts} \\
        ${groups_file} \\
        -o leafcutter_ds_${contrast_name}/leafcutter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leafcutter: \$(pip show leafcutter 2>/dev/null | grep Version | sed 's/Version: //' || echo 'unknown')
    END_VERSIONS
    """
}
