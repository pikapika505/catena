#!/usr/bin/env python3
"""
bin/annotate_celltypist.py
Automated cell type annotation using CellTypist.
Adds predicted_labels and majority_voting columns to adata.obs.
"""

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--h5ad",      required=True)
    parser.add_argument("--model",     default="Immune_All_Low.pkl",
                        help="CellTypist model name from model zoo")
    parser.add_argument("--outdir",    default=".")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print(f"[CellTypist] Sample: {args.sample_id} | Model: {args.model}")

    # ── Load adata ───────────────────────────────────────────
    adata = sc.read_h5ad(args.h5ad)
    print(f"[CellTypist] Loaded {adata.n_obs:,} cells")

    # CellTypist requires log1p-normalized counts (sum to ~1e4)
    # If raw is stored, use it
    if adata.raw is not None:
        adata_ct = adata.raw.to_adata().copy()
        sc.pp.normalize_total(adata_ct, target_sum=1e4)
        sc.pp.log1p(adata_ct)
    else:
        adata_ct = adata.copy()

    # ── Download + run model ─────────────────────────────────
    print(f"[CellTypist] Downloading model: {args.model}")
    models.download_models(model=args.model, force_update=False)
    model = models.Model.load(model=args.model)

    print("[CellTypist] Running prediction with majority voting...")
    predictions = celltypist.annotate(
        adata_ct,
        model         = model,
        majority_voting = True,
        over_clustering = "leiden" if "leiden" in adata.obs.columns else None
    )

    # Transfer predictions back to original adata
    adata.obs["celltypist_predicted"]       = predictions.predicted_labels.predicted_labels
    adata.obs["celltypist_majority_voting"] = predictions.predicted_labels.majority_voting
    adata.obs["celltypist_conf_score"]      = predictions.probability_matrix.max(axis=1).values

    n_types = adata.obs["celltypist_majority_voting"].nunique()
    print(f"[CellTypist] {n_types} cell types identified (majority voting)")

    # ── Save probability matrix ──────────────────────────────
    prob_df = predictions.probability_matrix.copy()
    prob_df.index = adata.obs_names
    prob_df["sample_id"]          = args.sample_id
    prob_df["predicted_label"]    = adata.obs["celltypist_predicted"].values
    prob_df["majority_voting"]    = adata.obs["celltypist_majority_voting"].values
    prob_df.to_csv(
        os.path.join(args.outdir, f"{args.sample_id}_celltypist_probs.csv")
    )

    # ── Plots ────────────────────────────────────────────────
    pdf_path = os.path.join(args.outdir, f"{args.sample_id}_celltypist_plots.pdf")
    with PdfPages(pdf_path) as pdf:

        # UMAP by majority voting
        if "X_umap" in adata.obsm:
            fig = sc.pl.umap(
                adata, color=["celltypist_majority_voting", "celltypist_conf_score"],
                title=[f"{args.sample_id} — CellTypist (majority voting)",
                       f"{args.sample_id} — Confidence score"],
                frameon=False, show=False, return_fig=True
            )
            pdf.savefig(fig); plt.close(fig)

        # Cell type composition bar chart
        comp = (adata.obs["celltypist_majority_voting"]
                .value_counts(normalize=True)
                .sort_values(ascending=True))
        fig, ax = plt.subplots(figsize=(8, max(4, len(comp)*0.35)))
        comp.plot.barh(ax=ax, color="steelblue")
        ax.set_xlabel("Fraction of cells")
        ax.set_title(f"{args.sample_id} — Cell type composition (CellTypist)")
        plt.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # Cluster → cell type assignment heatmap
        if "leiden" in adata.obs.columns:
            ct_matrix = pd.crosstab(
                adata.obs["leiden"],
                adata.obs["celltypist_majority_voting"],
                normalize="index"
            )
            fig, ax = plt.subplots(figsize=(max(8, ct_matrix.shape[1]*0.6),
                                            max(4, ct_matrix.shape[0]*0.5)))
            import seaborn as sns
            sns.heatmap(ct_matrix, cmap="Blues", ax=ax, linewidths=0.5,
                        annot=ct_matrix.shape[0] <= 20)
            ax.set_title(f"{args.sample_id} — Leiden cluster → CellTypist")
            ax.set_xlabel("Cell type")
            ax.set_ylabel("Leiden cluster")
            plt.tight_layout()
            pdf.savefig(fig); plt.close(fig)

    # ── Save annotated h5ad ──────────────────────────────────
    out_h5ad = os.path.join(args.outdir, f"{args.sample_id}_annotated.h5ad")
    adata.write_h5ad(out_h5ad)
    print(f"[CellTypist] Saved annotated adata → {out_h5ad}")
    print("[CellTypist] Done.")


if __name__ == "__main__":
    main()
