// modules/annotate_manual.nf
// ─────────────────────────────────────────────────────────────
// Marker-gene-based manual annotation applied to Seurat object
// Reads a user-supplied JSON of { "cell_type": ["gene1","gene2"] }
// Scores each cluster with AddModuleScore → picks top-scoring label
// ─────────────────────────────────────────────────────────────

process ANNOTATE_MANUAL {
    tag "${sample_id}"
    label 'medium'
    publishDir "${params.outdir}/${sample_id}/annotation/manual", mode: 'copy'

    conda "${moduleDir}/../envs/seurat.yaml"

    input:
    tuple val(sample_id), path(rds)

    output:
    tuple val(sample_id), path("${sample_id}_annotated.rds"), emit: rds
    path "${sample_id}_manual_scores.csv",                     emit: scores
    path "${sample_id}_manual_plots.pdf",                      emit: plots

    script:
    def markers_arg = params.marker_genes_file ?
        "--markers_file ${params.marker_genes_file}" : ""
    """
    Rscript ${projectDir}/bin/annotate_manual.R \\
        --sample_id  ${sample_id} \\
        --rds        ${rds} \\
        ${markers_arg} \\
        --outdir     .
    """
}
