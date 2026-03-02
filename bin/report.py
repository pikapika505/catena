#!/usr/bin/env python3
"""
bin/report.py
Aggregate per-sample benchmark metrics into a single HTML report.
"""

import argparse
import os
import glob
import pandas as pd
from datetime import datetime


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>scRNA/Spatial Pipeline Report</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: 'Segoe UI', sans-serif; background: #f5f7fa; color: #333; }}
  header {{ background: #1a237e; color: white; padding: 24px 40px; }}
  header h1 {{ font-size: 1.8em; font-weight: 300; }}
  header p  {{ font-size: 0.9em; opacity: 0.8; margin-top: 4px; }}
  .container {{ max-width: 1100px; margin: 40px auto; padding: 0 20px; }}
  .card {{ background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08);
           padding: 24px; margin-bottom: 28px; }}
  h2 {{ font-size: 1.2em; color: #1a237e; margin-bottom: 16px;
        border-bottom: 2px solid #e8eaf6; padding-bottom: 8px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; }}
  th {{ background: #e8eaf6; color: #1a237e; padding: 10px 14px;
        text-align: left; font-weight: 600; }}
  td {{ padding: 9px 14px; border-bottom: 1px solid #f0f0f0; }}
  tr:hover td {{ background: #f5f7fa; }}
  .badge {{ display: inline-block; padding: 3px 10px; border-radius: 12px;
            font-size: 0.8em; font-weight: 600; }}
  .good  {{ background: #e8f5e9; color: #2e7d32; }}
  .ok    {{ background: #fff8e1; color: #f57f17; }}
  .poor  {{ background: #ffebee; color: #c62828; }}
  .stat-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
               gap: 16px; }}
  .stat {{ background: #f5f7fa; border-radius: 6px; padding: 16px; text-align: center; }}
  .stat .value {{ font-size: 2em; font-weight: 700; color: #1a237e; }}
  .stat .label {{ font-size: 0.8em; color: #666; margin-top: 4px; }}
  footer {{ text-align: center; padding: 24px; color: #999; font-size: 0.85em; }}
</style>
</head>
<body>
<header>
  <h1>🧬 scRNA-seq / Spatial Transcriptomics Pipeline Report</h1>
  <p>Generated: {timestamp} &nbsp;|&nbsp; Samples: {n_samples}</p>
</header>
<div class="container">

  <div class="card">
    <h2>Run Summary</h2>
    <div class="stat-grid">
      <div class="stat"><div class="value">{n_samples}</div><div class="label">Samples processed</div></div>
      <div class="stat"><div class="value">{mean_ari}</div><div class="label">Mean ARI</div></div>
      <div class="stat"><div class="value">{mean_nmi}</div><div class="label">Mean NMI</div></div>
      <div class="stat"><div class="value">{mean_jaccard}</div><div class="label">Mean Marker Jaccard</div></div>
    </div>
  </div>

  <div class="card">
    <h2>Per-Sample Benchmarking</h2>
    <table>
      <thead>
        <tr>
          <th>Sample</th>
          <th>Cells</th>
          <th>Seurat clusters</th>
          <th>Scanpy clusters</th>
          <th>ARI</th>
          <th>NMI</th>
          <th>Marker Jaccard</th>
          <th>Verdict</th>
        </tr>
      </thead>
      <tbody>
        {rows}
      </tbody>
    </table>
  </div>

  <div class="card">
    <h2>Metric Interpretation</h2>
    <table>
      <thead><tr><th>Metric</th><th>Description</th><th>Good</th><th>Acceptable</th><th>Poor</th></tr></thead>
      <tbody>
        <tr><td><b>ARI</b></td><td>Adjusted Rand Index — cluster label agreement corrected for chance</td>
            <td>≥ 0.7</td><td>0.4–0.7</td><td>&lt; 0.4</td></tr>
        <tr><td><b>NMI</b></td><td>Normalized Mutual Information — shared information between clusterings</td>
            <td>≥ 0.7</td><td>0.4–0.7</td><td>&lt; 0.4</td></tr>
        <tr><td><b>Marker Jaccard</b></td><td>Mean overlap of top-50 marker genes per cluster</td>
            <td>≥ 0.3</td><td>0.15–0.3</td><td>&lt; 0.15</td></tr>
      </tbody>
    </table>
  </div>

</div>
<footer>scrna-spatial-pipeline &nbsp;|&nbsp; Seurat + Scanpy parallel analysis</footer>
</body>
</html>
"""


def score_badge(val, thresholds):
    """Return HTML badge based on value vs (good, ok) thresholds."""
    try:
        v = float(val)
    except (ValueError, TypeError):
        return f'<span class="badge ok">N/A</span>'
    good, ok = thresholds
    if v >= good:
        cls = "good"
    elif v >= ok:
        cls = "ok"
    else:
        cls = "poor"
    return f'<span class="badge {cls}">{v:.3f}</span>'


def verdict(ari, nmi):
    try:
        a, n = float(ari), float(nmi)
        if a >= 0.7 and n >= 0.7:
            return '<span class="badge good">High agreement ✓</span>'
        elif a >= 0.4 and n >= 0.4:
            return '<span class="badge ok">Moderate ⚠</span>'
        else:
            return '<span class="badge poor">Low agreement ✗</span>'
    except:
        return '<span class="badge ok">N/A</span>'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics_dir", default=".")
    parser.add_argument("--outdir",      default=".")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── Load all metrics CSVs ─────────────────────────────────
    csv_files = glob.glob(os.path.join(args.metrics_dir, "*_benchmark_metrics.csv"))
    if not csv_files:
        print("[Report] No benchmark metrics found — generating empty report")
        dfs = [pd.DataFrame(columns=[
            "sample_id","n_shared_barcodes","seurat_n_clusters",
            "scanpy_n_clusters","ARI","NMI","mean_marker_jaccard"])]
    else:
        dfs = [pd.read_csv(f) for f in csv_files]

    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(os.path.join(args.outdir, "summary_metrics.csv"), index=False)

    # ── Build table rows ──────────────────────────────────────
    rows = ""
    for _, row in df.iterrows():
        rows += f"""
        <tr>
          <td><b>{row.get('sample_id','')}</b></td>
          <td>{int(row.get('n_shared_barcodes',0)):,}</td>
          <td>{row.get('seurat_n_clusters','')}</td>
          <td>{row.get('scanpy_n_clusters','')}</td>
          <td>{score_badge(row.get('ARI','NA'), (0.7, 0.4))}</td>
          <td>{score_badge(row.get('NMI','NA'), (0.7, 0.4))}</td>
          <td>{score_badge(row.get('mean_marker_jaccard','NA'), (0.3, 0.15))}</td>
          <td>{verdict(row.get('ARI','NA'), row.get('NMI','NA'))}</td>
        </tr>"""

    # ── Summary stats ─────────────────────────────────────────
    def safe_mean(col):
        try:
            return f"{pd.to_numeric(df[col], errors='coerce').mean():.3f}"
        except:
            return "N/A"

    html = HTML_TEMPLATE.format(
        timestamp    = datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_samples    = len(df),
        mean_ari     = safe_mean("ARI"),
        mean_nmi     = safe_mean("NMI"),
        mean_jaccard = safe_mean("mean_marker_jaccard"),
        rows         = rows,
    )

    out_path = os.path.join(args.outdir, "pipeline_report.html")
    with open(out_path, "w") as f:
        f.write(html)

    print(f"[Report] HTML report → {out_path}")
    print("[Report] Done.")


if __name__ == "__main__":
    main()
