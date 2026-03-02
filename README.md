# scrna-spatial-pipeline

Parallel **Seurat** (R) + **Scanpy** (Python) preprocessing, clustering, annotation, and benchmarking pipeline for scRNA-seq and Visium spatial transcriptomics data.  
Built with **Nextflow DSL2** and **conda** environments — designed to run on the **Hellbender HPC** (SLURM) or locally.

```
CellRanger / SpaceRanger outs/
          │
          ▼
    [QC & Filtering]          ← adaptive (MAD) or fixed thresholds, or both
          │
    ┌─────┴─────┐
    ▼           ▼
[Seurat]    [Scanpy]          ← parallel branches, identical input
    │           │
    ▼           ▼
[Cluster]  [Leiden]
    │           │
    ▼           ▼
[Manual]  [CellTypist]        ← annotation (optional)
    │           │
    └─────┬─────┘
          ▼
    [Benchmarking]            ← ARI, NMI, marker Jaccard, cluster stats
          │
          ▼
    [HTML Report]
```

---

## Quick Start (Hellbender HPC)

### 1. One-time setup

```bash
# Load Nextflow (or install via conda)
module load nextflow   # or: conda install -c bioconda nextflow

# Clone the repo
git clone https://github.com/YOUR_USERNAME/scrna-spatial-pipeline.git
cd scrna-spatial-pipeline

# Point conda cache to scratch to avoid quota issues
echo 'conda.cacheDir = "/scratch/$USER/.nextflow/conda_cache"' >> ~/.nextflow/config
```

### 2. Edit your SLURM account

Open `nextflow.config` and set your SLURM account:
```groovy
process.clusterOptions = '--account=YOUR_SLURM_ACCOUNT'
```

### 3. Run — scRNA-seq

```bash
# Edit input path in params file
vim params/scrna_default.yaml

nextflow run main.nf \
  -profile hellbender \
  -params-file params/scrna_default.yaml \
  -resume
```

### 4. Run — Visium spatial

```bash
vim params/spatial_default.yaml

nextflow run main.nf \
  -profile hellbender \
  -params-file params/spatial_default.yaml \
  -resume
```

### 5. Run locally (for testing)

```bash
nextflow run main.nf \
  -profile local \
  --input_dir /path/to/cellranger_outputs \
  --outdir ./test_results \
  --mode scrna
```

---

## Input Structure

The pipeline expects one subdirectory per sample, following standard CellRanger/SpaceRanger output layout:

```
input_dir/
├── sample_A/
│   └── outs/
│       ├── filtered_feature_bc_matrix/   ← scRNA-seq
│       └── ...
├── sample_B/
│   └── outs/
│       ├── filtered_feature_bc_matrix/
│       └── spatial/                       ← Visium (spatial mode)
└── sample_C/
    └── outs/
```

---

## Key Parameters

All parameters can be set in a YAML file (`-params-file`) or passed on the command line (`--param value`).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mode` | `scrna` | `scrna` or `spatial` |
| `input_dir` | — | Path to directory containing sample subdirs |
| `outdir` | `./results` | Output directory |
| `qc_mode` | `adaptive` | `fixed`, `adaptive`, or `both` |
| `min_genes` | `200` | Fixed QC: min genes per cell |
| `max_genes` | `6000` | Fixed QC: max genes per cell |
| `max_pct_mt` | `20` | Fixed QC: max % mitochondrial |
| `mad_nmads` | `3` | Adaptive QC: number of MADs from median |
| `n_pcs` | `30` | Number of PCs for neighbors |
| `resolution` | `0.5` | Clustering resolution |
| `n_neighbors` | `20` | k-nearest neighbors |
| `annotate` | `true` | Run annotation step |
| `annotation_method` | `both` | `celltypist`, `manual`, or `both` |
| `celltypist_model` | `Immune_All_Low.pkl` | CellTypist model name |
| `marker_genes_file` | — | Path to JSON for manual annotation |

---

## Outputs

```
results/
├── pipeline_report.html          ← main summary report (open in browser)
├── summary_metrics.csv           ← aggregated benchmarking metrics
├── pipeline_info/
│   ├── report_*.html             ← Nextflow execution report
│   ├── timeline_*.html           ← job timeline
│   └── trace_*.tsv               ← resource usage per task
└── {sample_id}/
    ├── qc/
    │   ├── {sample}_qc_metrics.csv
    │   ├── {sample}_qc_plots.pdf
    │   └── {sample}_filtered/       ← filtered MTX + h5ad
    ├── seurat/
    │   ├── {sample}_seurat.rds
    │   ├── {sample}_seurat_meta.csv
    │   ├── {sample}_seurat_markers.csv
    │   └── {sample}_seurat_plots.pdf
    ├── scanpy/
    │   ├── {sample}_scanpy.h5ad
    │   ├── {sample}_scanpy_meta.csv
    │   ├── {sample}_scanpy_markers.csv
    │   └── {sample}_scanpy_plots.pdf
    ├── annotation/
    │   ├── celltypist/
    │   └── manual/
    └── benchmark/
        ├── {sample}_benchmark_metrics.csv
        └── {sample}_benchmark_plots.pdf
```

---

## Benchmarking Metrics

| Metric | Description | Interpretation |
|--------|-------------|----------------|
| **ARI** | Adjusted Rand Index | ≥ 0.7: high agreement ✓ |
| **NMI** | Normalized Mutual Information | ≥ 0.7: high agreement ✓ |
| **Marker Jaccard** | Mean overlap of top-50 marker genes per cluster | ≥ 0.3: good overlap |
| **Cluster size CV** | Coefficient of variation of cluster sizes | Lower = more balanced |

A high ARI/NMI with low Jaccard can indicate the tools find the same cell groupings but call different marker genes — worth investigating.

---

## Manual Annotation Marker Genes

Edit `params/example_markers.json` or supply your own JSON:

```json
{
  "T_cells":    ["CD3D", "CD3E", "CD8A"],
  "B_cells":    ["CD19", "MS4A1", "CD79A"],
  "Monocytes":  ["CD14", "LYZ", "S100A8"]
}
```

Pass it with: `--marker_genes_file params/your_markers.json`

---

## CellTypist Models

Available models: https://www.celltypist.org/models

Common choices:
- `Immune_All_Low.pkl` — immune cells, fine-grained
- `Immune_All_High.pkl` — immune cells, broad
- `Human_Lung_Atlas.pkl` — lung
- `Human_Colon_Cancer.pkl` — colon cancer
- `Developing_Human_Brain.pkl` — brain

---

## Resuming Failed Runs

Nextflow caches completed tasks. Add `-resume` to restart from where it stopped:

```bash
nextflow run main.nf -profile hellbender -params-file params/scrna_default.yaml -resume
```

---

## Dependencies

All software is managed automatically via conda. No manual installation needed.

- **Nextflow** ≥ 23.04 (`conda install -c bioconda nextflow`)
- **mamba** (optional, faster conda solves: `conda install -n base -c conda-forge mamba`)
- Python packages: `scanpy`, `squidpy`, `celltypist`, `scikit-learn`, `seaborn`
- R packages: `Seurat` ≥ 5.0, `optparse`, `patchwork`, `jsonlite`

---

## Citation / Acknowledgements

If you use this pipeline, please cite the underlying tools:
- **Seurat**: Hao et al., Cell 2021 / Hao et al., Nature Biotechnology 2024
- **Scanpy**: Wolf et al., Genome Biology 2018
- **Leiden**: Traag et al., Scientific Reports 2019
- **CellTypist**: Domínguez Conde et al., Science 2022
- **Nextflow**: Di Tommaso et al., Nature Biotechnology 2017
