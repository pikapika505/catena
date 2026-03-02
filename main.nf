#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
// scrna-spatial-pipeline
// Parallel Seurat + Scanpy preprocessing, clustering,
// annotation, and benchmarking of scRNA-seq / Visium data
// ============================================================

log.info """
╔══════════════════════════════════════════════════════════╗
║          scrna-spatial-pipeline  v1.0.0                  ║
║  Seurat ↔ Scanpy parallel analysis + benchmarking        ║
╚══════════════════════════════════════════════════════════╝
  mode          : ${params.mode}
  input_dir     : ${params.input_dir}
  outdir        : ${params.outdir}
  qc_mode       : ${params.qc_mode}
  annotate      : ${params.annotate}
  annotation_method : ${params.annotation_method}
  celltypist_model  : ${params.celltypist_model}
""".stripIndent()

// ── Import modules ───────────────────────────────────────────
include { QC             } from './modules/qc'
include { SEURAT         } from './modules/seurat'
include { SCANPY         } from './modules/scanpy'
include { ANNOTATE_CELLTYPIST } from './modules/annotate_celltypist'
include { ANNOTATE_MANUAL     } from './modules/annotate_manual'
include { BENCHMARK      } from './modules/benchmark'
include { REPORT         } from './modules/report'

// ── Validate inputs ──────────────────────────────────────────
def validateParams() {
    assert params.mode in ['scrna', 'spatial'] :
        "ERROR: --mode must be 'scrna' or 'spatial'. Got: ${params.mode}"
    assert params.qc_mode in ['fixed', 'adaptive', 'both'] :
        "ERROR: --qc_mode must be 'fixed', 'adaptive', or 'both'"
    assert params.annotation_method in ['celltypist', 'manual', 'both'] :
        "ERROR: --annotation_method must be 'celltypist', 'manual', or 'both'"
    def inputDir = file(params.input_dir)
    assert inputDir.exists() :
        "ERROR: input_dir does not exist: ${params.input_dir}"
}

// ── Main workflow ────────────────────────────────────────────
workflow {

    validateParams()

    // Collect sample directories (each subdir = one sample)
    // Expects: input_dir/sample_name/outs/  (CellRanger/SpaceRanger layout)
    ch_samples = Channel
        .fromPath("${params.input_dir}/*/outs", type: 'dir')
        .map { outs_dir ->
            def sample_id = outs_dir.parent.name
            tuple(sample_id, outs_dir)
        }

    ch_samples.view { sid, d -> "  → Found sample: ${sid}  [${d}]" }

    // ── Step 1: QC & filtering ───────────────────────────────
    QC(ch_samples)

    // ── Step 2: Parallel Seurat + Scanpy ────────────────────
    SEURAT(QC.out.filtered)
    SCANPY(QC.out.filtered)

    // ── Step 3: Annotation (optional) ───────────────────────
    if (params.annotate) {

        if (params.annotation_method in ['celltypist', 'both']) {
            ANNOTATE_CELLTYPIST(SCANPY.out.adata)
            ch_scanpy_final = ANNOTATE_CELLTYPIST.out.adata
        } else {
            ch_scanpy_final = SCANPY.out.adata
        }

        if (params.annotation_method in ['manual', 'both']) {
            ANNOTATE_MANUAL(SEURAT.out.rds)
            ch_seurat_final = ANNOTATE_MANUAL.out.rds
        } else {
            ch_seurat_final = SEURAT.out.rds
        }

    } else {
        ch_seurat_final = SEURAT.out.rds
        ch_scanpy_final = SCANPY.out.adata
    }

    // ── Step 4: Benchmarking ─────────────────────────────────
    // Join seurat + scanpy outputs by sample_id
    ch_bench_input = ch_seurat_final.join(ch_scanpy_final, by: 0)
    BENCHMARK(ch_bench_input)

    // ── Step 5: HTML Report ──────────────────────────────────
    REPORT(BENCHMARK.out.metrics.collect())
}

workflow.onComplete {
    log.info """
    ═══════════════════════════════════════════════
    Pipeline complete!
    Status   : ${workflow.success ? 'SUCCESS ✓' : 'FAILED ✗'}
    Duration : ${workflow.duration}
    Results  : ${params.outdir}
    ═══════════════════════════════════════════════
    """.stripIndent()
}
