// modules/benchmark.nf
// ─────────────────────────────────────────────────────────────
// Compare Seurat vs Scanpy clustering solutions
// Metrics: ARI, NMI, ASW (silhouette), marker gene overlap,
//          cluster size distributions, runtime/memory
// ─────────────────────────────────────────────────────────────

process BENCHMARK {
    tag "${sample_id}"
    label 'medium'
    publishDir "${params.outdir}/${sample_id}/benchmark", mode: 'copy'

    conda "${moduleDir}/../envs/scanpy.yaml"

    input:
    tuple val(sample_id), path(seurat_meta), path(scanpy_meta)

    output:
    tuple val(sample_id), path("${sample_id}_benchmark_metrics.csv"), emit: metrics
    path "${sample_id}_benchmark_plots.pdf",                           emit: plots

    script:
    """
    python ${projectDir}/bin/benchmark.py \\
        --sample_id   ${sample_id} \\
        --seurat_meta ${seurat_meta} \\
        --scanpy_meta ${scanpy_meta} \\
        --outdir      .
    """
}
