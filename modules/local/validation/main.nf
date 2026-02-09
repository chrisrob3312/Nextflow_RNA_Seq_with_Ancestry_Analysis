/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validation Module - QC checks and artifact detection
    Ensures data quality and flags potential artifacts across analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VALIDATE_DE_CONCORDANCE {
    label 'process_medium'

    input:
    path deseq2_results
    path limma_results
    val  contrast_name

    output:
    path "concordance_${contrast_name}/",         emit: results
    path "concordance_${contrast_name}/concordance_summary.tsv", emit: summary
    path "concordance_findings.md",               emit: findings
    path "versions.yml",                          emit: versions

    script:
    """
    mkdir -p concordance_${contrast_name}

    Rscript ${projectDir}/bin/validate_de_concordance.R \\
        --deseq2 ${deseq2_results} \\
        --limma ${limma_results} \\
        --contrast ${contrast_name} \\
        --output-dir concordance_${contrast_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

process VALIDATE_SEX_CHECK {
    label 'process_low'

    input:
    path count_matrix
    path metadata

    output:
    path "sex_check/",            emit: results
    path "sex_check_findings.md", emit: findings
    path "versions.yml",          emit: versions

    script:
    """
    mkdir -p sex_check

    Rscript ${projectDir}/bin/validate_sex_check.R \\
        --counts ${count_matrix} \\
        --metadata ${metadata} \\
        --output-dir sex_check

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

process VALIDATE_EXPRESSION_OUTLIERS {
    label 'process_medium'

    input:
    path normalized_counts
    path metadata

    output:
    path "expression_outliers/",            emit: results
    path "expression_outlier_findings.md",  emit: findings
    path "versions.yml",                    emit: versions

    script:
    """
    mkdir -p expression_outliers

    Rscript ${projectDir}/bin/validate_expression_outliers.R \\
        --counts ${normalized_counts} \\
        --metadata ${metadata} \\
        --output-dir expression_outliers \\
        --sd-threshold ${params.outlier_sd_threshold ?: 3} \\
        --min-correlation ${params.min_sample_correlation ?: 0.8}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

process VALIDATE_GENOMIC_INFLATION {
    label 'process_low'

    input:
    path de_results    // DE results with p-values
    val  contrast_name

    output:
    path "inflation_${contrast_name}/",     emit: results
    path "inflation_findings.md",           emit: findings
    path "versions.yml",                    emit: versions

    script:
    """
    mkdir -p inflation_${contrast_name}

    Rscript ${projectDir}/bin/validate_genomic_inflation.R \\
        --de-results ${de_results} \\
        --contrast ${contrast_name} \\
        --output-dir inflation_${contrast_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

process COMPILE_VALIDATION_REPORT {
    label 'process_low'

    input:
    path validation_results   // All validation output directories
    path findings_files       // All findings markdown files

    output:
    path "validation_report.md",   emit: report
    path "validation_summary.tsv", emit: summary
    path "versions.yml",           emit: versions

    script:
    """
    cat <<-HEADER > validation_report.md
    # Pipeline Validation Report
    **Generated:** \$(date -u '+%Y-%m-%d %H:%M:%S UTC')

    ## Overview
    This report summarizes quality control checks and validation steps
    performed during the analysis to identify potential artifacts.

    ---

    HEADER

    # Append all findings
    for f in ${findings_files}; do
        if [ -f "\$f" ]; then
            cat "\$f" >> validation_report.md
            echo "" >> validation_report.md
            echo "---" >> validation_report.md
            echo "" >> validation_report.md
        fi
    done

    # Create summary table
    echo -e "check\\tstatus\\tnotes" > validation_summary.tsv
    echo -e "de_concordance\\tSee report\\tDESeq2 vs limma-voom agreement" >> validation_summary.tsv
    echo -e "sex_check\\tSee report\\tXIST/Y-gene expression vs reported sex" >> validation_summary.tsv
    echo -e "expression_outliers\\tSee report\\tPCA/correlation-based outlier detection" >> validation_summary.tsv
    echo -e "genomic_inflation\\tSee report\\tLambda genomic inflation factor" >> validation_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1)
    END_VERSIONS
    """
}
