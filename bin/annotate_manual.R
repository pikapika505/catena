#!/usr/bin/env Rscript
# bin/annotate_manual.R
# Score each Seurat cluster against user-supplied marker gene lists.
# Marker genes JSON format: { "T cells": ["CD3D","CD3E","CD3G"], "B cells": [...] }
# Uses AddModuleScore → assigns each cluster the cell type with highest mean score.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(jsonlite)
  library(optparse)
})

# ── Args ─────────────────────────────────────────────────────
option_list <- list(
  make_option("--sample_id",    type="character"),
  make_option("--rds",          type="character"),
  make_option("--markers_file", type="character", default=NULL,
              help="JSON: { 'CellType': ['gene1','gene2', ...] }"),
  make_option("--outdir",       type="character", default=".")
)
opt <- parse_args(OptionParser(option_list=option_list))

message(sprintf("[Manual] Sample: %s", opt$sample_id))

# ── Load Seurat object ───────────────────────────────────────
obj <- readRDS(opt$rds)
message(sprintf("[Manual] Loaded %d cells, %d clusters",
                ncol(obj), nlevels(obj$seurat_clusters)))

# ── Load or define marker genes ─────────────────────────────
if (!is.null(opt$markers_file) && file.exists(opt$markers_file)) {
  marker_list <- fromJSON(opt$markers_file)
  message(sprintf("[Manual] Loaded %d cell type signatures from %s",
                  length(marker_list), opt$markers_file))
} else {
  # Fallback: a minimal general-purpose signature set
  message("[Manual] No markers_file provided — using built-in generic signatures")
  message("[Manual] TIP: Pass --markers_file your_markers.json for tissue-specific results")
  marker_list <- list(
    "T_cells"        = c("CD3D","CD3E","CD8A","CD4","TRAC"),
    "B_cells"        = c("CD19","MS4A1","CD79A","CD79B"),
    "NK_cells"       = c("NCAM1","NKG7","GNLY","KLRD1"),
    "Monocytes"      = c("CD14","LYZ","S100A8","S100A9","FCGR3A"),
    "Dendritic"      = c("FCER1A","CST3","IL3RA","CLEC4C"),
    "Plasma_cells"   = c("IGHG1","MZB1","SDC1","CD38"),
    "Epithelial"     = c("EPCAM","KRT18","KRT8","CDH1"),
    "Endothelial"    = c("PECAM1","VWF","CDH5","CLDN5"),
    "Fibroblasts"    = c("COL1A1","COL1A2","DCN","PDGFRA"),
    "Macrophages"    = c("CD68","MRC1","MARCO","ADGRE1")
  )
}

# ── AddModuleScore per cell type ─────────────────────────────
DefaultAssay(obj) <- ifelse("RNA" %in% names(obj@assays), "RNA", names(obj@assays)[1])

score_cols <- c()
for (ct in names(marker_list)) {
  genes <- intersect(marker_list[[ct]], rownames(obj))
  if (length(genes) == 0) {
    message(sprintf("  [skip] %s — no genes found in object", ct))
    next
  }
  col_name <- paste0("score_", gsub("[^A-Za-z0-9]", "_", ct))
  obj <- AddModuleScore(obj, features = list(genes),
                        name = col_name, ctrl = 50)
  # AddModuleScore appends "1" to the name
  actual_col <- paste0(col_name, "1")
  colnames(obj@meta.data)[colnames(obj@meta.data) == actual_col] <- col_name
  score_cols <- c(score_cols, col_name)
  message(sprintf("  [scored] %s (%d/%d genes present)",
                  ct, length(genes), length(marker_list[[ct]])))
}

# ── Assign cell type per cluster (max mean score) ────────────
if (length(score_cols) > 0) {
  score_mat <- obj@meta.data[, score_cols, drop=FALSE]
  colnames(score_mat) <- gsub("^score_", "", colnames(score_mat))

  cluster_scores <- score_mat %>%
    mutate(cluster = obj$seurat_clusters) %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean), .groups="drop")

  cluster_annot <- cluster_scores %>%
    tidyr::pivot_longer(-cluster, names_to="cell_type", values_to="score") %>%
    group_by(cluster) %>%
    slice_max(score, n=1) %>%
    ungroup() %>%
    select(cluster, predicted_celltype=cell_type, max_score=score)

  obj$manual_celltype <- cluster_annot$predicted_celltype[
    match(obj$seurat_clusters, cluster_annot$cluster)
  ]

  n_types <- length(unique(obj$manual_celltype))
  message(sprintf("[Manual] %d cell types assigned across %d clusters",
                  n_types, nlevels(obj$seurat_clusters)))
} else {
  warning("[Manual] No markers scored — manual_celltype set to 'Unknown'")
  obj$manual_celltype <- "Unknown"
  cluster_scores <- data.frame()
}

# ── Plots ────────────────────────────────────────────────────
pdf(sprintf("%s/%s_manual_plots.pdf", opt$outdir, opt$sample_id),
    width=12, height=8)

# UMAP colored by manual annotation
p1 <- DimPlot(obj, reduction="umap", group.by="manual_celltype",
              label=TRUE, repel=TRUE) +
  ggtitle(sprintf("%s — Manual annotation", opt$sample_id)) +
  theme_minimal()
print(p1)

# Side-by-side: clusters vs annotation
p2 <- DimPlot(obj, reduction="umap", group.by="seurat_clusters",
              label=TRUE, repel=TRUE) +
  ggtitle("Seurat clusters") + theme_minimal()
p3 <- DimPlot(obj, reduction="umap", group.by="manual_celltype",
              label=TRUE, repel=TRUE) +
  ggtitle("Manual annotation") + theme_minimal()
print(p2 | p3)

# Module score violin plots (top 6 cell types max)
if (length(score_cols) > 0) {
  plot_cols <- score_cols[seq_len(min(6, length(score_cols)))]
  p_vln <- VlnPlot(obj, features=plot_cols, group.by="seurat_clusters",
                   ncol=3, pt.size=0) &
    theme(axis.text.x=element_text(size=7))
  print(p_vln)
}

# Cluster score heatmap
if (nrow(cluster_scores) > 0) {
  score_m <- as.matrix(cluster_scores[,-1])
  rownames(score_m) <- paste0("C", cluster_scores$cluster)
  heatmap(score_m, scale="column",
          main=sprintf("%s — Cluster × Cell type scores", opt$sample_id),
          cexRow=0.9, cexCol=0.7)
}

dev.off()

# ── Export ───────────────────────────────────────────────────
scores_out <- obj@meta.data[, c("seurat_clusters","manual_celltype",
                                 score_cols), drop=FALSE]
scores_out$barcode    <- rownames(scores_out)
scores_out$sample_id  <- opt$sample_id

write.csv(scores_out,
          file=sprintf("%s/%s_manual_scores.csv", opt$outdir, opt$sample_id),
          row.names=FALSE)

rds_out <- sprintf("%s/%s_annotated.rds", opt$outdir, opt$sample_id)
saveRDS(obj, rds_out)
message(sprintf("[Manual] Saved annotated RDS → %s", rds_out))
message("[Manual] Done.")
