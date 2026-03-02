// modules/scanpy.nf
// ─────────────────────────────────────────────────────────────
// Scanpy branch: normalization → HVG → PCA → neighbors →
// Leiden clustering → UMAP → marker genes
// ─────────────────────────────────────────────────────────────

process SCANPY {
    tag "${sample_id}"
    label 'high'
    publishDir "${params.outdir}/${sample_id}/scanpy", mode: 'copy'

    conda "${moduleDir}/../envs/scanpy.yaml"

    input:
    tuple val(sample_id), path(filtered_dir)

    output:
    tuple val(sample_id), path("${sample_id}_scanpy.h5ad"),    emit: adata
    tuple val(sample_id), path("${sample_id}_scanpy_meta.csv"), emit: metadata
    path "${sample_id}_scanpy_markers.csv",                     emit: markers
    path "${sample_id}_scanpy_plots.pdf",                       emit: plots

    script:
    """
    python ${projectDir}/bin/scanpy_pipeline.py \\
        --sample_id    ${sample_id} \\
        --input_dir    ${filtered_dir} \\
        --mode         ${params.mode} \\
        --n_pcs        ${params.n_pcs} \\
        --resolution   ${params.resolution} \\
        --n_neighbors  ${params.n_neighbors} \\
        --outdir       .
    """
}
