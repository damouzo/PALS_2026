# Output

A successful run produces the following directory structure under
`${params.outdir}` (default: `results/`):

```
results/
├── preprocess_atac/
│   ├── all_peaks.csv             # chr-prefixed peaks, one per line
│   ├── cicero_connections.csv    # Peak1, Peak2, coaccess
│   ├── rna.h5ad                  # SeuratDisk export of the RNA assay
│   └── atac.h5ad                 # SeuratDisk export of the ATAC assay
│
├── infer_grn/
│   ├── base_GRN.parquet          # CellOracle base GRN (TF -> target)
│   ├── oracle.celloracle.oracle  # Serialized Oracle object
│   ├── links.celloracle.links    # Serialized Links object
│   ├── network_scores.xlsx       # Per-cluster network metrics
│   └── plots/
│       ├── vector_field_KO_GATA1.png   # ⭐ in silico perturbation
│       ├── vector_field_KO_MYB.png
│       ├── vector_field_KO_KLF1.png
│       ├── network_ranked_score.png
│       └── score_comparison.png
│
└── pipeline_info/
    ├── execution_trace.txt
    ├── execution_timeline.html
    ├── execution_report.html
    └── pipeline_dag.svg
```

## Highlights

- `infer_grn/plots/vector_field_KO_<TF>.png` — the visual centerpiece of the
  pipeline. Each plot shows a UMAP with a vector field that predicts how cells
  would shift if the named TF were knocked out.
- `infer_grn/base_GRN.parquet` — the inferred base Gene Regulatory Network. The
  primary input for any downstream TF analysis.
- `pipeline_info/execution_timeline.html` — Nextflow's automatic task-duration
  timeline. Useful for resource accounting and reports.
