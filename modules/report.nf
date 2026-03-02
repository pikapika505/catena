// modules/report.nf
// ─────────────────────────────────────────────────────────────
// Collects all per-sample metrics and generates a single
// summary HTML report for the whole run
// ─────────────────────────────────────────────────────────────

process REPORT {
    label 'low'
    publishDir "${params.outdir}", mode: 'copy'

    conda "${moduleDir}/../envs/scanpy.yaml"

    input:
    path(metrics_files)  // list of *_benchmark_metrics.csv

    output:
    path "pipeline_report.html", emit: html
    path "summary_metrics.csv",  emit: summary

    script:
    """
    python ${projectDir}/bin/report.py \\
        --metrics_dir . \\
        --outdir      .
    """
}
