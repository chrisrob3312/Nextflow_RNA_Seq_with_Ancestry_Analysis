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
    path herv_results, stageAs: 'herv_results/*'

    output:
    path "merged_neoantigens.tsv",      emit: merged
    path "neoantigen_summary.tsv",      emit: summary
    path "neoantigen_burden.tsv",       emit: burden
    path "neoantigen_findings.md",      emit: findings
    path "versions.yml",                emit: versions

    script:
    def herv_arg = herv_results.name != 'NO_HERV' ? "--herv-dir herv_results" : ''
    """
    python3 ${projectDir}/bin/neoantigen_filter.py \\
        --pvacseq-dir . \\
        --neofuse-dir . \\
        --snaf-dir . \\
        ${herv_arg} \\
        --binding-threshold ${params.binding_threshold} \\
        --output-merged merged_neoantigens.tsv \\
        --output-summary neoantigen_summary.tsv \\
        --output-burden neoantigen_burden.tsv

    # Findings
    TOTAL=\$(awk 'NR>1' merged_neoantigens.tsv 2>/dev/null | wc -l || echo 0)
    cat <<-FINDINGS > neoantigen_findings.md
    ### Merged Neoantigen Summary
    - Total unique neoantigen candidates: \${TOTAL}
    - Sources: pVACseq (SNV/indel), NeoFuse (fusion), SNAF (splicing)${herv_arg ? ', HERV/TE' : ''}
    - See neoantigen_burden.tsv for per-sample burden
    FINDINGS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process TELESCOPE_HERV {
    tag "$meta.sample_id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path herv_annotation

    output:
    tuple val(meta), path("*-telescope_report.tsv"), emit: herv_counts
    path "*.findings.md", emit: findings

    script:
    """
    telescope assign \\
        --attribute locus \\
        --no_feature_key __no_feature \\
        --ncpu ${task.cpus} \\
        --outdir . \\
        ${bam} \\
        ${herv_annotation}

    mv telescope-telescope_report.tsv ${meta.sample_id}-telescope_report.tsv

    n_expressed=\$(awk -F'\\t' 'NR>1 && \$3 > 0' ${meta.sample_id}-telescope_report.tsv | wc -l)
    cat <<-FINDINGS > ${meta.sample_id}.herv.findings.md
    ## HERV/TE Expression - ${meta.sample_id}
    - HERVs/TEs with non-zero expression: \${n_expressed}
    FINDINGS
    """
}

process HERVQUANT {
    tag "$meta.sample_id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path hervquant_ref

    output:
    tuple val(meta), path("*.hervquant.tsv"), emit: hervquant_results

    script:
    """
    # hervQuant quantifies HERV expression from RNA-seq
    hervQuant.py \\
        --bam ${bam} \\
        --ref ${hervquant_ref} \\
        --output ${meta.sample_id}.hervquant.tsv \\
        --threads ${task.cpus}
    """
}

process NEOANTIGEN_BURDEN_ANALYSIS {
    label 'process_medium'

    input:
    path merged_neoantigens
    path metadata
    path ancestry_proportions

    output:
    path "neoantigen_burden/", emit: burden_results
    path "*.findings.md", emit: findings

    script:
    """
    mkdir -p neoantigen_burden
    Rscript <<'REOF'
    library(data.table)
    library(ggplot2)

    neo <- fread("${merged_neoantigens}")
    meta <- fread("${metadata}")
    anc <- tryCatch(fread("${ancestry_proportions}"), error = function(e) NULL)

    if (!is.null(anc)) meta <- merge(meta, anc, by = "sample_id", all.x = TRUE)

    # Per-sample neoantigen burden by source
    burden <- neo[, .(
        total_neoantigens = .N,
        snv_neoantigens = sum(source == "snv", na.rm = TRUE),
        fusion_neoantigens = sum(source == "fusion", na.rm = TRUE),
        splicing_neoantigens = sum(source == "splicing", na.rm = TRUE)
    ), by = sample_id]

    burden <- merge(burden, meta, by = "sample_id", all.x = TRUE)

    # --- Relapse vs non-relapse comparison ---
    if ("relapse_status" %in% names(burden) && length(unique(na.omit(burden\$relapse_status))) >= 2) {
        relapse_test <- wilcox.test(total_neoantigens ~ relapse_status, data = burden)
        sink("neoantigen_burden/relapse_comparison.txt")
        cat("Neoantigen burden: Relapse vs Non-relapse\\n")
        cat(sprintf("Relapse median: %.1f\\n", median(burden[relapse_status == "relapse"]\$total_neoantigens, na.rm = TRUE)))
        cat(sprintf("Non-relapse median: %.1f\\n", median(burden[relapse_status == "no_relapse"]\$total_neoantigens, na.rm = TRUE)))
        cat(sprintf("Wilcoxon p-value: %.4e\\n", relapse_test\$p.value))
        sink()

        p1 <- ggplot(burden[!is.na(relapse_status)], aes(x = relapse_status, y = total_neoantigens, fill = relapse_status)) +
            geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.5) +
            labs(title = "Neoantigen Burden by Relapse Status", x = "", y = "Total Neoantigens") +
            theme_bw() + theme(legend.position = "none")
        ggsave("neoantigen_burden/burden_by_relapse.pdf", p1, width = 6, height = 5)
    }

    # --- MRD status comparison ---
    if ("mrd_status" %in% names(burden) && length(unique(na.omit(burden\$mrd_status))) >= 2) {
        mrd_test <- wilcox.test(total_neoantigens ~ mrd_status, data = burden)
        sink("neoantigen_burden/mrd_comparison.txt")
        cat("Neoantigen burden: MRD-positive vs MRD-negative\\n")
        cat(sprintf("MRD-positive median: %.1f\\n", median(burden[mrd_status == "positive"]\$total_neoantigens, na.rm = TRUE)))
        cat(sprintf("MRD-negative median: %.1f\\n", median(burden[mrd_status == "negative"]\$total_neoantigens, na.rm = TRUE)))
        cat(sprintf("Wilcoxon p-value: %.4e\\n", mrd_test\$p.value))
        sink()

        p2 <- ggplot(burden[!is.na(mrd_status)], aes(x = mrd_status, y = total_neoantigens, fill = mrd_status)) +
            geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.5) +
            labs(title = "Neoantigen Burden by MRD Status", x = "", y = "Total Neoantigens") +
            theme_bw() + theme(legend.position = "none")
        ggsave("neoantigen_burden/burden_by_mrd.pdf", p2, width = 6, height = 5)
    }

    # --- Ancestry group comparison ---
    if ("graf_category" %in% names(burden) && length(unique(na.omit(burden\$graf_category))) >= 2) {
        anc_test <- kruskal.test(total_neoantigens ~ graf_category, data = burden[!is.na(graf_category)])
        sink("neoantigen_burden/ancestry_comparison.txt")
        cat("Neoantigen burden by ancestry group (Kruskal-Wallis)\\n")
        for (grp in unique(na.omit(burden\$graf_category))) {
            cat(sprintf("%s median: %.1f (n=%d)\\n", grp,
                median(burden[graf_category == grp]\$total_neoantigens, na.rm = TRUE),
                sum(burden\$graf_category == grp, na.rm = TRUE)))
        }
        cat(sprintf("Kruskal-Wallis p-value: %.4e\\n", anc_test\$p.value))
        sink()

        p3 <- ggplot(burden[!is.na(graf_category)], aes(x = graf_category, y = total_neoantigens, fill = graf_category)) +
            geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.5) +
            labs(title = "Neoantigen Burden by Ancestry", x = "", y = "Total Neoantigens") +
            theme_bw() + theme(legend.position = "none")
        ggsave("neoantigen_burden/burden_by_ancestry.pdf", p3, width = 8, height = 5)
    }

    # --- Per-source breakdown stacked bar ---
    burden_long <- melt(burden, id.vars = "sample_id",
                        measure.vars = c("snv_neoantigens", "fusion_neoantigens", "splicing_neoantigens"),
                        variable.name = "source", value.name = "count")
    burden_long[, source := gsub("_neoantigens", "", source)]

    p4 <- ggplot(burden_long, aes(x = reorder(sample_id, -count), y = count, fill = source)) +
        geom_bar(stat = "identity") +
        labs(title = "Neoantigen Burden by Source", x = "Sample", y = "Neoantigen Count") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
    ggsave("neoantigen_burden/burden_by_source.pdf", p4, width = max(8, nrow(burden) * 0.3), height = 5)

    # Save burden table
    fwrite(burden, "neoantigen_burden/neoantigen_burden_annotated.tsv", sep = "\\t")
    REOF

    cat <<-FINDINGS > neoantigen_burden.findings.md
    ## Neoantigen Burden Analysis
    - Samples analyzed: \$(tail -n +2 neoantigen_burden/neoantigen_burden_annotated.tsv | wc -l)
    - Stratified by: relapse status, MRD status, ancestry
    FINDINGS
    """
}

process ANTIGEN_DB_CROSSREF {
    label 'process_low'

    input:
    path merged_neoantigens
    path antigen_db

    output:
    path "antigen_crossref/", emit: crossref_results

    script:
    """
    mkdir -p antigen_crossref
    python3 <<'PYEOF'
import pandas as pd
import os

neo = pd.read_csv("${merged_neoantigens}", sep="\\t")
db_path = "${antigen_db}"

known_antigens = set()
if os.path.exists(db_path):
    db = pd.read_csv(db_path, sep="\\t")
    if "antigen" in db.columns:
        known_antigens = set(db["antigen"].str.upper())
    elif "gene" in db.columns:
        known_antigens = set(db["gene"].str.upper())

if "gene" in neo.columns and known_antigens:
    neo["known_cancer_antigen"] = neo["gene"].str.upper().isin(known_antigens)
    known_hits = neo[neo["known_cancer_antigen"]]
    known_hits.to_csv("antigen_crossref/known_antigen_neoantigens.tsv", sep="\\t", index=False)

    with open("antigen_crossref/crossref_summary.txt", "w") as f:
        f.write(f"Total neoantigens: {len(neo)}\\n")
        f.write(f"In known cancer antigen DB: {len(known_hits)}\\n")
        f.write(f"Known antigen genes: {', '.join(sorted(known_hits['gene'].unique()))}\\n")
else:
    with open("antigen_crossref/crossref_summary.txt", "w") as f:
        f.write("No antigen database provided or no gene column in neoantigens\\n")

neo.to_csv("antigen_crossref/neoantigens_annotated.tsv", sep="\\t", index=False)
PYEOF
    """
}
