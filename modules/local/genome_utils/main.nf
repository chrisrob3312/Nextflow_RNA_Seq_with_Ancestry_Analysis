/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Genome Utilities Module
    - Auto-detect genome build from BAM header (hg19/GRCh37 vs hg38/GRCh38)
    - CrossMap liftover for BAM files between builds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DETECT_GENOME_BUILD {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path(bam), path(bai), env(BUILD), emit: bam_with_build
    path "${meta.id}.genome_build.txt",                 emit: build_info
    path "versions.yml",                                emit: versions

    script:
    """
    #!/bin/bash
    set -euo pipefail

    BUILD="unknown"

    # Extract header
    samtools view -H ${bam} > header.txt

    # Check for chromosome naming conventions
    HAS_CHR=\$(grep -c '^@SQ.*SN:chr' header.txt || true)
    HAS_NOCHR=\$(grep -c '^@SQ.*SN:[0-9]' header.txt || true)

    # Get chr1/1 length
    # hg38/GRCh38: chr1 = 248956422
    # hg19/GRCh37: chr1 = 249250621
    CHR1_LEN=\$(grep -P '^@SQ\\tSN:(chr)?1\\t' header.txt | grep -oP 'LN:\\K[0-9]+' || echo "0")

    if [ "\${CHR1_LEN}" -eq 248956422 ]; then
        BUILD="hg38"
    elif [ "\${CHR1_LEN}" -eq 249250621 ]; then
        if [ "\${HAS_CHR}" -gt 0 ]; then
            BUILD="hg19"
        else
            BUILD="GRCh37"
        fi
    else
        # Fallback: check for known assembly identifiers in header
        if grep -q 'GRCh38\\|hg38\\|GCA_000001405.15' header.txt; then
            BUILD="hg38"
        elif grep -q 'GRCh37\\|hg19\\|GCA_000001405.1[^5]' header.txt; then
            BUILD="hg19"
        fi
    fi

    # Chromosome style
    if [ "\${HAS_CHR}" -gt 0 ]; then
        CHR_STYLE="UCSC"
    elif [ "\${HAS_NOCHR}" -gt 0 ]; then
        CHR_STYLE="Ensembl"
    else
        CHR_STYLE="unknown"
    fi

    echo -e "sample_id\\tgenome_build\\tchr_style\\tchr1_length" > ${meta.id}.genome_build.txt
    echo -e "${meta.id}\\t\${BUILD}\\t\${CHR_STYLE}\\t\${CHR1_LEN}" >> ${meta.id}.genome_build.txt

    # BUILD env variable is captured by Nextflow via env() output qualifier

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    BUILD="hg38"
    echo -e "sample_id\\tgenome_build\\tchr_style\\tchr1_length" > ${meta.id}.genome_build.txt
    echo -e "${meta.id}\\thg38\\tUCSC\\t248956422" >> ${meta.id}.genome_build.txt
    touch versions.yml
    """
}

process CROSSMAP_BAM {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path  target_fasta      // Target genome FASTA (e.g., hg38)
    path  chain_file        // Chain file (e.g., hg19ToHg38.over.chain.gz)

    output:
    tuple val(meta), path("${meta.id}.liftover.sorted.bam"), path("${meta.id}.liftover.sorted.bam.bai"), emit: bam
    path "${meta.id}.liftover.unmap.bam",  emit: unmapped, optional: true
    path "${meta.id}.liftover_stats.txt",  emit: stats
    path "versions.yml",                   emit: versions

    script:
    """
    # CrossMap BAM liftover - handles spliced alignments correctly
    CrossMap.py bam \\
        ${chain_file} \\
        ${bam} \\
        ${meta.id}.liftover.bam \\
        -a ${target_fasta}

    # Sort the lifted-over BAM
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${meta.id}.liftover.sorted.bam \\
        ${meta.id}.liftover.bam

    samtools index ${meta.id}.liftover.sorted.bam

    # Collect liftover statistics
    TOTAL=\$(samtools view -c ${bam})
    MAPPED=\$(samtools view -c -F 4 ${meta.id}.liftover.sorted.bam)
    UNMAPPED=\$(samtools view -c ${meta.id}.liftover.unmap.bam 2>/dev/null || echo "0")
    PCT_MAPPED=\$(echo "scale=2; \${MAPPED} * 100 / \${TOTAL}" | bc)

    echo -e "sample_id\\ttotal_reads\\tmapped_after_liftover\\tunmapped\\tpct_mapped" > ${meta.id}.liftover_stats.txt
    echo -e "${meta.id}\\t\${TOTAL}\\t\${MAPPED}\\t\${UNMAPPED}\\t\${PCT_MAPPED}" >> ${meta.id}.liftover_stats.txt

    rm -f ${meta.id}.liftover.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crossmap: \$(CrossMap.py 2>&1 | grep -i version | sed 's/.*version //' || echo 'unknown')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}

process VALIDATE_GENOME_BUILDS {
    label 'process_low'

    input:
    path build_files  // All .genome_build.txt files

    output:
    path "genome_build_summary.tsv", emit: summary
    path "versions.yml",             emit: versions

    script:
    """
    # Combine all build detection results
    head -1 \$(ls *.genome_build.txt | head -1) > genome_build_summary.tsv
    for f in *.genome_build.txt; do
        tail -n +2 "\$f" >> genome_build_summary.tsv
    done

    # Check if all samples have the same build
    BUILDS=\$(tail -n +2 genome_build_summary.tsv | cut -f2 | sort -u)
    N_BUILDS=\$(echo "\$BUILDS" | wc -l)

    if [ "\$N_BUILDS" -gt 1 ]; then
        echo "WARNING: Mixed genome builds detected: \$BUILDS" >&2
        echo "Liftover will be required for consistent analysis." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1)
    END_VERSIONS
    """
}
