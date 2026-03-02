#!/usr/bin/env python3
"""
bin/benchmark.py
Compare Seurat vs Scanpy clustering solutions.
Metrics: ARI, NMI, ASW (silhouette), marker gene Jaccard overlap,
         cluster count, cluster size CV, runtime.
"""

import argparse
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.metrics.pairwise import cosine_similarity


def jaccard(set_a: set, set_b: set) -> float:
    if not set_a and not set_b:
        return 1.0
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)


def marker_overlap(seurat_markers_path: str, scanpy_markers_path: str,
                   top_n: int = 50) -> dict:
    """Per-cluster Jaccard of top-N marker genes."""
    srt = pd.read_csv(seurat_markers_path)
    scp = pd.read_csv(scanpy_markers_path)

    # Seurat: columns cluster, gene, avg_log2FC
    # Scanpy: columns group, names, logfoldchanges
    srt_top = (srt.sort_values("avg_log2FC", ascending=False)
                  .groupby("cluster")
                  .head(top_n)
                  .groupby("cluster")["gene"]
                  .apply(set)
                  .to_dict())

    scp_top = (scp.sort_values("logfoldchanges", ascending=False)
                  .groupby("group")
                  .head(top_n)
                  .groupby("group")["names"]
                  .apply(set)
                  .to_dict())

    # Match by rank (cluster 0 ↔ cluster 0 etc.)
    all_clusters = sorted(set(srt_top.keys()) | {str(k) for k in scp_top.keys()})
    scores = {}
    for c in all_clusters:
        s_genes = srt_top.get(c, set())
        q_genes = scp_top.get(str(c), scp_top.get(int(c) if str(c).isdigit() else c, set()))
        scores[str(c)] = round(jaccard(s_genes, q_genes), 4)

    return scores


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id",   required=True)
    parser.add_argument("--seurat_meta", required=True)
    parser.add_argument("--scanpy_meta", required=True)
    parser.add_argument("--outdir",      default=".")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print(f"[Benchmark] Sample: {args.sample_id}")

    # ── Load metadata ─────────────────────────────────────────
    srt = pd.read_csv(args.seurat_meta).set_index("barcode")
    scp = pd.read_csv(args.scanpy_meta).set_index("barcode")

    # Align on shared barcodes
    shared = srt.index.intersection(scp.index)
    n_shared = len(shared)
    print(f"[Benchmark] Shared barcodes: {n_shared:,} "
          f"(Seurat: {len(srt):,}, Scanpy: {len(scp):,})")

    srt = srt.loc[shared]
    scp = scp.loc[shared]

    srt_clusters = srt["tool_cluster"].astype(str)
    scp_clusters = scp["tool_cluster"].astype(str)

    # ── Clustering agreement metrics ──────────────────────────
    ari = adjusted_rand_score(srt_clusters, scp_clusters)
    nmi = normalized_mutual_info_score(srt_clusters, scp_clusters,
                                       average_method="arithmetic")
    n_srt = srt_clusters.nunique()
    n_scp = scp_clusters.nunique()

    print(f"[Benchmark] ARI  : {ari:.4f}")
    print(f"[Benchmark] NMI  : {nmi:.4f}")
    print(f"[Benchmark] Clusters — Seurat: {n_srt} | Scanpy: {n_scp}")

    # ── Cluster size distribution stats ──────────────────────
    def cluster_cv(clusters: pd.Series) -> float:
        sizes = clusters.value_counts()
        return sizes.std() / sizes.mean() if sizes.mean() > 0 else 0.0

    cv_srt = cluster_cv(srt_clusters)
    cv_scp = cluster_cv(scp_clusters)

    # ── Marker gene Jaccard (if files exist) ─────────────────
    srt_markers_path = args.seurat_meta.replace("_meta.csv", "_markers.csv")
    scp_markers_path = args.scanpy_meta.replace("_meta.csv", "_markers.csv")
    marker_scores = {}
    mean_jaccard  = np.nan

    if os.path.exists(srt_markers_path) and os.path.exists(scp_markers_path):
        try:
            marker_scores = marker_overlap(srt_markers_path, scp_markers_path)
            mean_jaccard  = np.mean(list(marker_scores.values()))
            print(f"[Benchmark] Mean marker Jaccard: {mean_jaccard:.4f}")
        except Exception as e:
            print(f"[Benchmark] WARNING: marker overlap failed: {e}")

    # ── Compose metrics table ─────────────────────────────────
    metrics = {
        "sample_id":         args.sample_id,
        "n_shared_barcodes": n_shared,
        "seurat_n_clusters": n_srt,
        "scanpy_n_clusters": n_scp,
        "ARI":               round(ari, 4),
        "NMI":               round(nmi, 4),
        "seurat_cluster_size_cv":  round(cv_srt, 4),
        "scanpy_cluster_size_cv":  round(cv_scp, 4),
        "mean_marker_jaccard":     round(mean_jaccard, 4) if not np.isnan(mean_jaccard) else "NA",
    }

    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(
        os.path.join(args.outdir, f"{args.sample_id}_benchmark_metrics.csv"),
        index=False
    )

    # ── Plots ─────────────────────────────────────────────────
    pdf_path = os.path.join(args.outdir, f"{args.sample_id}_benchmark_plots.pdf")
    with PdfPages(pdf_path) as pdf:

        # 1. Summary bar chart of key metrics
        fig, ax = plt.subplots(figsize=(7, 4))
        metric_vals = {"ARI": ari, "NMI": nmi}
        if not np.isnan(mean_jaccard):
            metric_vals["Marker Jaccard"] = mean_jaccard
        bars = ax.bar(list(metric_vals.keys()), list(metric_vals.values()),
                      color=["#2196F3","#4CAF50","#FF9800"][:len(metric_vals)])
        ax.set_ylim(0, 1.15)
        ax.set_ylabel("Score")
        ax.set_title(f"{args.sample_id} — Seurat vs Scanpy agreement")
        for bar, val in zip(bars, metric_vals.values()):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                    f"{val:.3f}", ha="center", fontsize=11, fontweight="bold")
        ax.axhline(0.7, color="gray", linestyle="--", linewidth=0.8,
                   label="0.7 threshold")
        ax.legend(fontsize=9)
        plt.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # 2. Cluster size distributions
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        for ax, (label, clusters) in zip(axes, [
                ("Seurat", srt_clusters), ("Scanpy", scp_clusters)]):
            sizes = clusters.value_counts().sort_index()
            ax.bar(sizes.index.astype(str), sizes.values,
                   color="#2196F3" if label=="Seurat" else "#4CAF50")
            ax.set_xlabel("Cluster")
            ax.set_ylabel("Number of cells")
            ax.set_title(f"{label} cluster sizes  (CV={cluster_cv(clusters):.2f})")
            ax.tick_params(axis="x", rotation=45)
        fig.suptitle(f"{args.sample_id} — Cluster size distributions", fontsize=13)
        plt.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # 3. Confusion matrix (Seurat rows, Scanpy cols)
        conf = pd.crosstab(srt_clusters, scp_clusters, normalize="index")
        fig, ax = plt.subplots(figsize=(max(6, conf.shape[1]*0.5),
                                        max(5, conf.shape[0]*0.5)))
        sns.heatmap(conf, cmap="Blues", ax=ax, linewidths=0.3,
                    fmt=".2f",
                    annot=conf.shape[0] <= 20 and conf.shape[1] <= 20)
        ax.set_xlabel("Scanpy cluster")
        ax.set_ylabel("Seurat cluster")
        ax.set_title(f"{args.sample_id} — Cluster assignment overlap\n"
                     f"ARI={ari:.3f}  NMI={nmi:.3f}")
        plt.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # 4. Per-cluster marker Jaccard (if available)
        if marker_scores:
            clusters = sorted(marker_scores.keys(),
                               key=lambda x: int(x) if x.isdigit() else x)
            scores   = [marker_scores[c] for c in clusters]
            fig, ax  = plt.subplots(figsize=(max(6, len(clusters)*0.5), 4))
            ax.bar(clusters, scores, color="#9C27B0")
            ax.set_xlabel("Cluster")
            ax.set_ylabel("Jaccard index (top-50 markers)")
            ax.set_ylim(0, 1)
            ax.axhline(mean_jaccard, color="red", linestyle="--",
                       label=f"Mean={mean_jaccard:.3f}")
            ax.legend()
            ax.set_title(f"{args.sample_id} — Per-cluster marker gene overlap")
            plt.tight_layout()
            pdf.savefig(fig); plt.close(fig)

    print(f"[Benchmark] Metrics saved → {args.outdir}")
    print("[Benchmark] Done.")


if __name__ == "__main__":
    main()
