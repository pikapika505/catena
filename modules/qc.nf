// modules/qc.nf
// ─────────────────────────────────────────────────────────────
// QC and filtering — runs once per sample before splitting
// into Seurat / Scanpy branches. Outputs a filtered 10x MTX
// directory that both tools can read natively.
// ─────────────────────────────────────────────────────────────

process QC {
    tag "${sample_id}"
    label 'medium'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    conda "${moduleDir}/../envs/scanpy.yaml"

    input:
    tuple val(sample_id), path(outs_dir)

    output:
    tuple val(sample_id), path("${sample_id}_filtered"), emit: filtered
    path "${sample_id}_qc_metrics.csv",                  emit: metrics
    path "${sample_id}_qc_plots.pdf",                    emit: plots

    script:
    def adaptive_flag = params.qc_mode in ['adaptive', 'both'] ? '--adaptive' : ''
    def fixed_flag    = params.qc_mode in ['fixed', 'both']    ? '--fixed'    : ''
    """
    python ${projectDir}/bin/qc.py \\
        --sample_id       ${sample_id} \\
        --outs_dir        ${outs_dir} \\
        --mode            ${params.mode} \\
        --outdir          ${sample_id}_filtered \\
        ${adaptive_flag} \\
        ${fixed_flag} \\
        --min_genes       ${params.min_genes} \\
        --max_genes       ${params.max_genes} \\
        --max_pct_mt      ${params.max_pct_mt} \\
        --min_cells       ${params.min_cells} \\
        --mad_nmads       ${params.mad_nmads}
    """
}
