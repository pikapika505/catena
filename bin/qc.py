#!/usr/bin/env python3
"""
bin/qc.py
QC filtering for scRNA-seq or Visium spatial data.
Supports fixed thresholds, adaptive (MAD-based), or both.
Outputs a filtered MTX directory readable by both Seurat and Scanpy.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def mad(x: np.ndarray) -> float:
    return np.median(np.abs(x - np.median(x)))


def adaptive_thresholds(adata: sc.AnnData, nmads: int = 3) -> dict:
    """Compute per-metric MAD-based outlier thresholds."""
    thresholds = {}
    for metric in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
        vals = adata.obs[metric].values
        med  = np.median(vals)
        m    = mad(vals)
        thresholds[metric] = {
            "lower": med - nmads * m,
            "upper": med + nmads * m,
        }
    return thresholds


def flag_cells(adata: sc.AnnData, thresholds: dict) -> pd.Series:
    """Return boolean mask: True = keep."""
    keep = pd.Series(True, index=adata.obs_names)
    keep &= adata.obs["n_genes_by_counts"] >= thresholds["n_genes_by_counts"]["lower"]
    keep &= adata.obs["n_genes_by_counts"] <= thresholds["n_genes_by_counts"]["upper"]
    keep &= adata.obs["total_counts"]       <= thresholds["total_counts"]["upper"]
    keep &= adata.obs["pct_counts_mt"]      <= thresholds["pct_counts_mt"]["upper"]
    return keep


def plot_qc(adata_raw: sc.AnnData, adata_filt: sc.AnnData,
            thresholds: dict, sample_id: str, pdf: PdfPages, mode_label: str):
    """Violin + scatter plots before and after filtering."""
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(f"{sample_id} — QC  ({mode_label})", fontsize=14)

    for col, metric in enumerate(metrics):
        for row, (label, adata) in enumerate([("Before", adata_raw), ("After", adata_filt)]):
            ax = axes[row, col]
            vals = adata.obs[metric]
            ax.violinplot(vals, showmedians=True)
            ax.set_title(f"{label}: {metric.replace('_', ' ')}")
            ax.set_xticks([])

            # Draw threshold lines
            if metric in thresholds:
                lo = thresholds[metric].get("lower", None)
                hi = thresholds[metric].get("upper", None)
                if lo is not None and lo > vals.min():
                    ax.axhline(lo, color="red", linestyle="--", linewidth=0.8)
                if hi is not None and hi < vals.max():
                    ax.axhline(hi, color="red", linestyle="--", linewidth=0.8)

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # Scatter: total_counts vs n_genes, colored by pct_mt
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (label, adata) in zip(axes, [("Before", adata_raw), ("After", adata_filt)]):
        sc_plot = ax.scatter(
            adata.obs["total_counts"],
            adata.obs["n_genes_by_counts"],
            c=adata.obs["pct_counts_mt"],
            s=1, alpha=0.5, cmap="viridis"
        )
        plt.colorbar(sc_plot, ax=ax, label="% MT")
        ax.set_xlabel("Total counts")
        ax.set_ylabel("Genes detected")
        ax.set_title(f"{label} filtering  (n={adata.n_obs:,})")
    fig.suptitle(f"{sample_id} — Count vs Genes scatter", fontsize=13)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="QC filtering for scRNA-seq / Visium data")
    parser.add_argument("--sample_id",   required=True)
    parser.add_argument("--outs_dir",    required=True, help="CellRanger/SpaceRanger outs/ directory")
    parser.add_argument("--mode",        default="scrna", choices=["scrna", "spatial"])
    parser.add_argument("--outdir",      required=True)
    parser.add_argument("--adaptive",    action="store_true")
    parser.add_argument("--fixed",       action="store_true")
    parser.add_argument("--min_genes",   type=int,   default=200)
    parser.add_argument("--max_genes",   type=int,   default=6000)
    parser.add_argument("--max_pct_mt",  type=float, default=20.0)
    parser.add_argument("--min_cells",   type=int,   default=3)
    parser.add_argument("--mad_nmads",   type=int,   default=3)
    args = parser.parse_args()

    if not args.adaptive and not args.fixed:
        print("WARNING: Neither --adaptive nor --fixed specified. Defaulting to adaptive.")
        args.adaptive = True

    os.makedirs(args.outdir, exist_ok=True)
    sc.settings.verbosity = 1

    # ── Load data ────────────────────────────────────────────
    print(f"[QC] Loading {args.mode} data from {args.outs_dir}")
    if args.mode == "scrna":
        mtx_path = os.path.join(args.outs_dir, "filtered_feature_bc_matrix")
        adata = sc.read_10x_mtx(mtx_path, var_names="gene_symbols", cache=True)
    else:  # spatial
        adata = sq.read.visium(args.outs_dir)
        adata.var_names_make_unique()

    print(f"[QC] Loaded: {adata.n_obs:,} barcodes × {adata.n_vars:,} genes")

    # ── Basic QC metrics ─────────────────────────────────────
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    adata_raw = adata.copy()

    # ── Build thresholds ─────────────────────────────────────
    thresholds = {}
    mode_labels = []

    if args.fixed:
        thresholds["fixed"] = {
            "n_genes_by_counts": {"lower": args.min_genes,  "upper": args.max_genes},
            "total_counts":      {"lower": 0,               "upper": np.inf},
            "pct_counts_mt":     {"lower": 0,               "upper": args.max_pct_mt},
        }
        mode_labels.append("fixed")

    if args.adaptive:
        thresholds["adaptive"] = adaptive_thresholds(adata, nmads=args.mad_nmads)
        mode_labels.append(f"adaptive (MAD×{args.mad_nmads})")

    # When both are requested, intersect (take stricter threshold per metric)
    if args.fixed and args.adaptive:
        combined = {}
        for metric in thresholds["fixed"]:
            combined[metric] = {
                "lower": max(thresholds["fixed"][metric]["lower"],
                             thresholds["adaptive"][metric]["lower"]),
                "upper": min(thresholds["fixed"][metric]["upper"],
                             thresholds["adaptive"][metric]["upper"]),
            }
        final_thresholds = combined
        mode_label = "fixed + adaptive (intersection)"
    elif args.adaptive:
        final_thresholds = thresholds["adaptive"]
        mode_label = f"adaptive (MAD×{args.mad_nmads})"
    else:
        final_thresholds = thresholds["fixed"]
        mode_label = "fixed"

    print(f"[QC] Thresholds ({mode_label}):")
    for metric, bounds in final_thresholds.items():
        lo = f"{bounds['lower']:.1f}"
        hi = f"{bounds['upper']:.1f}" if bounds["upper"] != np.inf else "∞"
        print(f"     {metric}: [{lo}, {hi}]")

    # ── Filter ───────────────────────────────────────────────
    keep = flag_cells(adata, final_thresholds)
    adata = adata[keep].copy()
    n_removed = adata_raw.n_obs - adata.n_obs
    print(f"[QC] Retained {adata.n_obs:,} / {adata_raw.n_obs:,} barcodes "
          f"({n_removed:,} removed, {100*n_removed/adata_raw.n_obs:.1f}%)")

    # ── Save QC metrics CSV ──────────────────────────────────
    metrics_df = adata_raw.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
    metrics_df["passed_qc"] = keep.values
    metrics_df["sample_id"] = args.sample_id
    metrics_df.to_csv(f"{args.sample_id}_qc_metrics.csv")

    # ── Save QC plots ────────────────────────────────────────
    with PdfPages(f"{args.sample_id}_qc_plots.pdf") as pdf:
        plot_qc(adata_raw, adata, final_thresholds, args.sample_id, pdf, mode_label)

    # ── Write filtered MTX (Seurat + Scanpy compatible) ──────
    print(f"[QC] Writing filtered MTX to {args.outdir}")
    sc.pp.filter_genes(adata, min_cells=0)  # ensure matrix is clean
    adata.write_h5ad(os.path.join(args.outdir, "filtered.h5ad"))

    # Also write 10x MTX format for Seurat
    import scipy.io
    mtx_out = os.path.join(args.outdir, "matrix")
    os.makedirs(mtx_out, exist_ok=True)
    scipy.io.mmwrite(os.path.join(mtx_out, "matrix.mtx"), adata.X.T)
    adata.obs_names.to_frame().to_csv(
        os.path.join(mtx_out, "barcodes.tsv"), header=False, index=False)
    adata.var_names.to_frame().assign(
        feature_type="Gene Expression"
    ).to_csv(os.path.join(mtx_out, "features.tsv"),
             sep="\t", header=False, index=True)

    # Copy spatial metadata if spatial mode
    if args.mode == "spatial" and hasattr(adata, "uns") and "spatial" in adata.uns:
        import shutil
        spatial_src = os.path.join(args.outs_dir, "spatial")
        spatial_dst = os.path.join(args.outdir, "spatial")
        if os.path.exists(spatial_src) and not os.path.exists(spatial_dst):
            shutil.copytree(spatial_src, spatial_dst)

    print("[QC] Done.")


if __name__ == "__main__":
    main()
