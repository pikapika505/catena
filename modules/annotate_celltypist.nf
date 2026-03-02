// modules/annotate_celltypist.nf
// ─────────────────────────────────────────────────────────────
// Automated cell type annotation with CellTypist
// Takes Scanpy .h5ad, adds predicted_labels + majority_voting
// ─────────────────────────────────────────────────────────────

process ANNOTATE_CELLTYPIST {
    tag "${sample_id}"
    label 'medium'
    publishDir "${params.outdir}/${sample_id}/annotation/celltypist", mode: 'copy'

    conda "${moduleDir}/../envs/scanpy.yaml"

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), path("${sample_id}_annotated.h5ad"), emit: adata
    path "${sample_id}_celltypist_probs.csv",                   emit: probs
    path "${sample_id}_celltypist_plots.pdf",                   emit: plots

    script:
    """
    python ${projectDir}/bin/annotate_celltypist.py \\
        --sample_id  ${sample_id} \\
        --h5ad       ${h5ad} \\
        --model      ${params.celltypist_model} \\
        --outdir     .
    """
}
