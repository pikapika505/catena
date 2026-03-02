#!/usr/bin/env Rscript
# bin/seurat_pipeline.R
# Normalization → HVG → PCA → Neighbors → Clustering → UMAP → Markers
# Handles both scRNA-seq and Visium spatial modes

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(optparse)
})

# ── Args ─────────────────────────────────────────────────────
option_list <- list(
  make_option("--sample_id",   type="character"),
  make_option("--input_dir",   type="character", help="Filtered dir from qc.py"),
  make_option("--mode",        type="character", default="scrna"),
  make_option("--n_pcs",       type="integer",   default=30),
  make_option("--resolution",  type="double",    default=0.5),
  make_option("--n_neighbors", type="integer",   default=20),
  make_option("--outdir",      type="character", default=".")
)
opt <- parse_args(OptionParser(option_list=option_list))

message(sprintf("[Seurat] Sample: %s | Mode: %s", opt$sample_id, opt$mode))

# ── Load data ────────────────────────────────────────────────
mtx_path <- file.path(opt$input_dir, "matrix")

if (opt$mode == "scrna") {
  counts <- Read10X(data.dir = mtx_path)
  obj <- CreateSeuratObject(counts = counts, project = opt$sample_id)
  message(sprintf("[Seurat] Loaded %d cells × %d genes",
                  ncol(obj), nrow(obj)))

} else {  # spatial
  # SpaceRanger layout: parent of 'matrix' contains 'spatial/'
  spatial_dir <- file.path(opt$input_dir)
  obj <- Load10X_Spatial(
    data.dir   = spatial_dir,
    filename   = "matrix/matrix.mtx",
    assay      = "Spatial"
  )
  message(sprintf("[Seurat] Loaded %d spots × %d genes",
                  ncol(obj), nrow(obj)))
}

# ── Normalization ────────────────────────────────────────────
if (opt$mode == "scrna") {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 1e4, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                               nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
} else {
  obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
}

# ── Dimensionality reduction ─────────────────────────────────
obj <- RunPCA(obj, npcs = opt$n_pcs, verbose = FALSE)

# Elbow plot
pdf_file <- sprintf("%s/%s_seurat_plots.pdf", opt$outdir, opt$sample_id)
pdf(pdf_file, width=10, height=8)

print(ElbowPlot(obj, ndims = opt$n_pcs) +
  ggtitle(sprintf("%s — Seurat Elbow Plot", opt$sample_id)))

# ── Neighbors + Clustering ───────────────────────────────────
obj <- FindNeighbors(obj, dims = 1:opt$n_pcs,
                     k.param = opt$n_neighbors, verbose = FALSE)
obj <- FindClusters(obj, resolution = opt$resolution, verbose = FALSE)
message(sprintf("[Seurat] Found %d clusters (resolution=%.2f)",
                nlevels(obj$seurat_clusters), opt$resolution))

# ── UMAP ────────────────────────────────────────────────────
obj <- RunUMAP(obj, dims = 1:opt$n_pcs, verbose = FALSE)

# UMAP plot
p1 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE) +
  ggtitle(sprintf("%s — Seurat clusters (res=%.2f)", opt$sample_id, opt$resolution)) +
  theme_minimal()
print(p1)

# Spatial plot (if spatial mode)
if (opt$mode == "spatial" && !is.null(obj@images)) {
  p_spatial <- SpatialDimPlot(obj, label = TRUE, label.size = 3) +
    ggtitle(sprintf("%s — Spatial clusters", opt$sample_id))
  print(p_spatial)
}

# ── Marker genes ─────────────────────────────────────────────
message("[Seurat] Finding marker genes...")
obj <- PrepSCTFindMarkers(obj, verbose = FALSE) # no-op for LogNorm
markers <- FindAllMarkers(
  obj,
  only.pos      = TRUE,
  min.pct       = 0.25,
  logfc.threshold = 0.5,
  verbose       = FALSE
)

# Top 5 markers per cluster heatmap
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

if (nrow(top5) > 0) {
  p_heat <- DoHeatmap(obj, features = top5$gene) +
    ggtitle(sprintf("%s — Top marker genes", opt$sample_id)) +
    theme(axis.text.y = element_text(size = 6))
  print(p_heat)
}

dev.off()

# ── Export metadata ──────────────────────────────────────────
meta <- obj@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA",
                           "seurat_clusters"), drop=FALSE]
# Rename for benchmarking compatibility
meta <- meta %>%
  rename(tool_cluster = seurat_clusters) %>%
  mutate(tool = "seurat", sample_id = opt$sample_id)
meta$barcode <- rownames(meta)

write.csv(meta,
          file = sprintf("%s/%s_seurat_meta.csv", opt$outdir, opt$sample_id),
          row.names = FALSE)

# Export markers
write.csv(markers,
          file = sprintf("%s/%s_seurat_markers.csv", opt$outdir, opt$sample_id),
          row.names = FALSE)

# ── Save RDS ─────────────────────────────────────────────────
rds_path <- sprintf("%s/%s_seurat.rds", opt$outdir, opt$sample_id)
saveRDS(obj, rds_path)
message(sprintf("[Seurat] Saved RDS → %s", rds_path))
message("[Seurat] Done.")
