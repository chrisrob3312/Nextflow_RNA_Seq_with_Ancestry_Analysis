/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Differential Splicing Module
    Tools: rMATS, Leafcutter, SplAdder + Bisbee
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
    path "rmats_${contrast_name}/",    emit: results
    path "rmats_findings.md",          emit: findings
    path "versions.yml",               emit: versions

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

    # Generate findings markdown
    cat <<-FINDINGS > rmats_findings.md
    ### rMATS Differential Splicing: ${contrast_name}
    **Event counts (FDR < 0.05):**
    FINDINGS
    for event_type in SE A3SS A5SS MXE RI; do
        FILE="rmats_${contrast_name}/\${event_type}.MATS.JC.txt"
        if [ -f "\$FILE" ]; then
            SIG=\$(awk -F'\\t' 'NR>1 && \$20<0.05' "\$FILE" | wc -l)
            TOTAL=\$(awk 'NR>1' "\$FILE" | wc -l)
            echo "- \${event_type}: \${SIG} significant / \${TOTAL} tested" >> rmats_findings.md
        fi
    done

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
    path "leafcutter_findings.md",          emit: findings
    path "versions.yml",                    emit: versions

    script:
    """
    mkdir -p leafcutter_ds_${contrast_name}

    Rscript \$(which leafcutter_ds.R) \\
        --num_threads ${task.cpus} \\
        ${counts} \\
        ${groups_file} \\
        -o leafcutter_ds_${contrast_name}/leafcutter

    # Generate findings
    SIG_CLUSTERS=\$(awk -F'\\t' 'NR>1 && \$5<0.05' leafcutter_ds_${contrast_name}/leafcutter_cluster_significance.txt 2>/dev/null | wc -l || echo 0)
    cat <<-FINDINGS > leafcutter_findings.md
    ### Leafcutter Differential Splicing: ${contrast_name}
    - Significant clusters (FDR < 0.05): \${SIG_CLUSTERS}
    FINDINGS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leafcutter: \$(pip show leafcutter 2>/dev/null | grep Version | sed 's/Version: //' || echo 'unknown')
    END_VERSIONS
    """
}

// ========================================
// SplAdder: Splicing graph construction + event detection
// ========================================

process SPLADDER_BUILD {
    label 'process_high'

    input:
    path bam_files   // All BAM files
    path bai_files   // All BAI files
    path gtf

    output:
    path "spladder_out/",          emit: graph
    path "spladder_out/merge_graphs_*.pickle", emit: event_files
    path "versions.yml",           emit: versions

    script:
    def confidence = params.spladder_confidence ?: 3
    def event_types = params.spladder_event_types ?: 'exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip'
    """
    # Create BAM file list
    ls *.bam > bam_list.txt

    # Build splicing graphs and detect alternative splicing events
    # SplAdder identifies: exon skipping (ES), intron retention (IR),
    # alternative 3'/5' splice sites (A3/A5), mutually exclusive exons (MXE)
    spladder build \\
        --bams bam_list.txt \\
        --annotation ${gtf} \\
        --outdir spladder_out \\
        --confidence ${confidence} \\
        --event-types ${event_types} \\
        --parallel ${task.cpus} \\
        --merge-strat merge_graphs \\
        --quantify-graph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spladder: \$(spladder --version 2>&1 | head -1 | sed 's/.*v//' || echo 'unknown')
    END_VERSIONS
    """
}

process SPLADDER_TEST {
    label 'process_high'

    input:
    path spladder_dir        // SplAdder output directory
    path condition_file      // Two-column: sample_id, condition (for comparison)
    val  contrast_name

    output:
    path "spladder_test_${contrast_name}/",   emit: results
    path "spladder_findings.md",              emit: findings
    path "versions.yml",                      emit: versions

    script:
    def event_types = params.spladder_event_types ?: 'exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip'
    """
    mkdir -p spladder_test_${contrast_name}

    # Run SplAdder differential testing (negative binomial GLM)
    for event in \$(echo ${event_types} | tr ',' ' '); do
        spladder test \\
            --outdir ${spladder_dir} \\
            --conditionA \$(awk -F'\\t' '\$2=="group1" {print \$1}' ${condition_file} | tr '\\n' ',') \\
            --conditionB \$(awk -F'\\t' '\$2=="group2" {print \$1}' ${condition_file} | tr '\\n' ',') \\
            --event-type \${event} \\
            --test-result spladder_test_${contrast_name}/\${event}_results.tsv \\
            --no-cap-exp-outliers || true
    done

    # Generate findings summary
    cat <<-FINDINGS > spladder_findings.md
    ### SplAdder Differential Splicing: ${contrast_name}
    **Event counts (p < 0.05, NB GLM):**
    FINDINGS
    for event in \$(echo ${event_types} | tr ',' ' '); do
        FILE="spladder_test_${contrast_name}/\${event}_results.tsv"
        if [ -f "\$FILE" ]; then
            SIG=\$(awk -F'\\t' 'NR>1 && \$NF<0.05' "\$FILE" | wc -l)
            TOTAL=\$(awk 'NR>1' "\$FILE" | wc -l)
            echo "- \${event}: \${SIG} significant / \${TOTAL} tested" >> spladder_findings.md
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spladder: \$(spladder --version 2>&1 | head -1 | sed 's/.*v//' || echo 'unknown')
    END_VERSIONS
    """
}

// ========================================
// Bisbee: Beta-binomial differential splicing & protein effects
// Depends on SplAdder event/count output
// ========================================

process BISBEE_PREP {
    label 'process_medium'

    input:
    path spladder_dir    // SplAdder output directory with HDF5 count files

    output:
    path "bisbee_input/",    emit: prepped
    path "versions.yml",     emit: versions

    script:
    """
    mkdir -p bisbee_input

    # Prepare SplAdder output for Bisbee
    # Extracts junction inclusion/exclusion counts from SplAdder HDF5 files
    python3 ${projectDir}/bin/bisbee_prep.py \\
        --spladder-dir ${spladder_dir} \\
        --output-dir bisbee_input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process BISBEE_DIFF {
    label 'process_high'

    input:
    path bisbee_input    // Prepped Bisbee input directory
    path condition_file  // Two-column: sample_id, condition
    val  contrast_name

    output:
    path "bisbee_diff_${contrast_name}/",   emit: results
    path "bisbee_diff_findings.md",         emit: findings
    path "versions.yml",                    emit: versions

    script:
    """
    mkdir -p bisbee_diff_${contrast_name}

    # Bisbee differential splicing using beta-binomial model
    # Accounts for overdispersion in junction count data
    Rscript ${projectDir}/bin/bisbee_diff.R \\
        --input-dir ${bisbee_input} \\
        --conditions ${condition_file} \\
        --contrast ${contrast_name} \\
        --output-dir bisbee_diff_${contrast_name} \\
        --threads ${task.cpus}

    # Findings
    SIG=\$(awk -F'\\t' 'NR>1 && \$NF<0.05' bisbee_diff_${contrast_name}/differential_splicing.tsv 2>/dev/null | wc -l || echo 0)
    cat <<-FINDINGS > bisbee_diff_findings.md
    ### Bisbee Differential Splicing: ${contrast_name}
    - Significant events (FDR < 0.05, beta-binomial): \${SIG}
    - Model: Beta-binomial GLM accounting for overdispersion
    FINDINGS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

process BISBEE_PROT {
    label 'process_medium'

    input:
    path bisbee_diff_results  // Bisbee differential splicing results
    path gtf
    path fasta

    output:
    path "bisbee_protein_effects/",         emit: protein_effects
    path "bisbee_protein_findings.md",      emit: findings
    path "versions.yml",                    emit: versions

    script:
    """
    mkdir -p bisbee_protein_effects

    # Predict protein-level effects of differentially spliced events
    # Maps splice events to reading frame changes, NMD targets, domain disruptions
    python3 ${projectDir}/bin/bisbee_prot.py \\
        --diff-results ${bisbee_diff_results} \\
        --gtf ${gtf} \\
        --fasta ${fasta} \\
        --output-dir bisbee_protein_effects

    # Findings
    NMD=\$(awk -F'\\t' 'NR>1 && \$0~/NMD/' bisbee_protein_effects/protein_effects.tsv 2>/dev/null | wc -l || echo 0)
    FRAMESHIFT=\$(awk -F'\\t' 'NR>1 && \$0~/frameshift/' bisbee_protein_effects/protein_effects.tsv 2>/dev/null | wc -l || echo 0)
    cat <<-FINDINGS > bisbee_protein_findings.md
    ### Bisbee Protein Effect Predictions
    - NMD-targeted events: \${NMD}
    - Frameshift events: \${FRAMESHIFT}
    - See bisbee_protein_effects/ for full domain/motif disruption analysis
    FINDINGS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process BISBEE_OUTLIER {
    label 'process_medium'

    input:
    path bisbee_input    // Prepped Bisbee input

    output:
    path "bisbee_outliers/",            emit: outliers
    path "bisbee_outlier_findings.md",  emit: findings
    path "versions.yml",                emit: versions

    script:
    """
    mkdir -p bisbee_outliers

    # Detect sample-level splicing outliers using beta-binomial model
    # Identifies per-sample events with unusually high/low inclusion
    Rscript ${projectDir}/bin/bisbee_outlier.R \\
        --input-dir ${bisbee_input} \\
        --output-dir bisbee_outliers \\
        --fdr-threshold ${params.bisbee_outlier_fdr ?: 0.05}

    # Findings
    cat <<-FINDINGS > bisbee_outlier_findings.md
    ### Bisbee Splice Outlier Detection
    FINDINGS
    if [ -f "bisbee_outliers/outlier_summary.tsv" ]; then
        echo "- Outlier samples and event counts:" >> bisbee_outlier_findings.md
        awk -F'\\t' 'NR>1 {print "  - " \$1 ": " \$2 " outlier events"}' bisbee_outliers/outlier_summary.tsv >> bisbee_outlier_findings.md
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}
