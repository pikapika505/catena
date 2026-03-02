// modules/seurat.nf
// ─────────────────────────────────────────────────────────────
// Seurat branch: normalization → HVG → PCA → neighbors →
// clustering → UMAP → marker genes
// ─────────────────────────────────────────────────────────────

process SEURAT {
    tag "${sample_id}"
    label 'high'
    publishDir "${params.outdir}/${sample_id}/seurat", mode: 'copy'

    conda "${moduleDir}/../envs/seurat.yaml"

    input:
    tuple val(sample_id), path(filtered_dir)

    output:
    tuple val(sample_id), path("${sample_id}_seurat.rds"),   emit: rds
    tuple val(sample_id), path("${sample_id}_seurat_meta.csv"), emit: metadata
    path "${sample_id}_seurat_markers.csv",                    emit: markers
    path "${sample_id}_seurat_plots.pdf",                      emit: plots

    script:
    """
    Rscript ${projectDir}/bin/seurat_pipeline.R \\
        --sample_id    ${sample_id} \\
        --input_dir    ${filtered_dir} \\
        --mode         ${params.mode} \\
        --n_pcs        ${params.n_pcs} \\
        --resolution   ${params.resolution} \\
        --n_neighbors  ${params.n_neighbors} \\
        --outdir       .
    """
}
