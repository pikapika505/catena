#!/usr/bin/env python3
"""
bin/scanpy_pipeline.py
Normalization → HVG → PCA → Neighbors → Leiden clustering → UMAP → Markers
Handles both scRNA-seq and Visium spatial modes.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import squidpy as sq
    HAS_SQUIDPY = True
except ImportError:
    HAS_SQUIDPY = False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id",   required=True)
    parser.add_argument("--input_dir",   required=True)
    parser.add_argument("--mode",        default="scrna", choices=["scrna", "spatial"])
    parser.add_argument("--n_pcs",       type=int,   default=30)
    parser.add_argument("--resolution",  type=float, default=0.5)
    parser.add_argument("--n_neighbors", type=int,   default=20)
    parser.add_argument("--outdir",      default=".")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sc.settings.verbosity = 2
    sc.settings.figdir = args.outdir

    print(f"[Scanpy] Sample: {args.sample_id} | Mode: {args.mode}")

    # ── Load filtered data ────────────────────────────────────
    h5ad_path = os.path.join(args.input_dir, "filtered.h5ad")
    if os.path.exists(h5ad_path):
        adata = sc.read_h5ad(h5ad_path)
        print(f"[Scanpy] Loaded h5ad: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    else:
        mtx_path = os.path.join(args.input_dir, "matrix")
        if args.mode == "scrna":
            adata = sc.read_10x_mtx(mtx_path, var_names="gene_symbols", cache=True)
        else:
            if not HAS_SQUIDPY:
                raise ImportError("squidpy required for spatial mode")
            adata = sq.read.visium(args.input_dir)
        adata.var_names_make_unique()
        print(f"[Scanpy] Loaded MTX: {adata.n_obs:,} × {adata.n_vars:,}")

    # ── Normalization ─────────────────────────────────────────
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # freeze normalized counts before HVG scaling

    # ── HVG ──────────────────────────────────────────────────
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=False)
    n_hvg = adata.var["highly_variable"].sum()
    print(f"[Scanpy] {n_hvg} highly variable genes selected")

    # Scale only HVGs for PCA
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)

    # ── PCA ──────────────────────────────────────────────────
    sc.tl.pca(adata_hvg, n_comps=args.n_pcs, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"]    = adata_hvg.uns["pca"]

    # ── Neighbors + Leiden clustering ────────────────────────
    sc.pp.neighbors(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)
    sc.tl.leiden(adata, resolution=args.resolution, key_added="leiden")
    n_clusters = adata.obs["leiden"].nunique()
    print(f"[Scanpy] Found {n_clusters} Leiden clusters (resolution={args.resolution})")

    # ── UMAP ─────────────────────────────────────────────────
    sc.tl.umap(adata)

    # ── Marker genes (rank_genes_groups) ─────────────────────
    print("[Scanpy] Computing marker genes...")
    sc.tl.rank_genes_groups(
        adata, groupby="leiden",
        method="wilcoxon", use_raw=True,
        key_added="rank_genes"
    )
    markers_df = sc.get.rank_genes_groups_df(
        adata, group=None, key="rank_genes",
        pval_cutoff=0.05, log2fc_min=0.5
    )
    markers_df["sample_id"] = args.sample_id

    # ── Plots ────────────────────────────────────────────────
    pdf_path = os.path.join(args.outdir, f"{args.sample_id}_scanpy_plots.pdf")
    with PdfPages(pdf_path) as pdf:

        # Variance explained
        fig, ax = plt.subplots(figsize=(8, 4))
        var_ratio = adata.uns["pca"]["variance_ratio"]
        ax.plot(np.arange(1, len(var_ratio)+1), np.cumsum(var_ratio)*100, "o-")
        ax.set_xlabel("PC")
        ax.set_ylabel("Cumulative variance explained (%)")
        ax.set_title(f"{args.sample_id} — PCA variance (Scanpy)")
        ax.axvline(args.n_pcs, color="red", linestyle="--", alpha=0.5,
                   label=f"Used: {args.n_pcs} PCs")
        ax.legend()
        plt.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # UMAP by cluster
        fig = sc.pl.umap(
            adata, color=["leiden"],
            title=f"{args.sample_id} — Leiden (res={args.resolution})",
            legend_loc="on data", frameon=False, show=False,
            return_fig=True
        )
        pdf.savefig(fig); plt.close(fig)

        # Top marker dot plot
        top_markers = (
            markers_df.sort_values("logfoldchanges", ascending=False)
            .groupby("group")
            .head(3)["names"]
            .unique()
            .tolist()
        )
        if top_markers:
            fig = sc.pl.dotplot(
                adata, var_names=top_markers[:30], groupby="leiden",
                title=f"{args.sample_id} — Top markers",
                show=False, return_fig=True
            )
            pdf.savefig(fig); plt.close(fig)

        # Spatial plot
        if args.mode == "spatial" and HAS_SQUIDPY and "spatial" in adata.uns:
            fig, ax = plt.subplots(figsize=(6, 6))
            sq.pl.spatial_scatter(adata, color="leiden", ax=ax,
                                  title=f"{args.sample_id} — Spatial clusters")
            pdf.savefig(fig); plt.close(fig)

    # ── Export ───────────────────────────────────────────────
    meta = adata.obs[["leiden"]].copy()
    meta.columns = ["tool_cluster"]
    meta["tool"]      = "scanpy"
    meta["sample_id"] = args.sample_id
    meta.index.name   = "barcode"
    meta.reset_index().to_csv(
        os.path.join(args.outdir, f"{args.sample_id}_scanpy_meta.csv"),
        index=False
    )

    markers_df.to_csv(
        os.path.join(args.outdir, f"{args.sample_id}_scanpy_markers.csv"),
        index=False
    )

    h5ad_out = os.path.join(args.outdir, f"{args.sample_id}_scanpy.h5ad")
    adata.write_h5ad(h5ad_out)
    print(f"[Scanpy] Saved → {h5ad_out}")
    print("[Scanpy] Done.")


if __name__ == "__main__":
    main()
