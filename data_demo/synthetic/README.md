# `data_demo/synthetic/` — Ultra-tiny synthetic test object (CI smoke tests)

This directory holds an ultra-tiny Seurat v4 object (≈ 300 KB, 150–200 cells)
that is **only** used for CI smoke tests and developer sanity-checks where the
multi-gigabyte biological dataset (`data_demo/processed/microglia_dam_demo.rds`)
is not available (no internet, no Docker, etc.).

## Regeneration

```bash
# Generate the synthetic object
Rscript data_demo/synthetic/build_synthetic_mini_seurat.R

# Sanity-check it
Rscript data_demo/synthetic/validate_synthetic_mini_seurat.R
```

The output is written to `data_demo/synthetic/mini_seurat.rds`. The script
tries several data sources in priority order:

1. `SeuratData::pbmcMultiome` (real 10x multiome, ~60 cells per CellType).
2. CellOracle example data (PBMC).
3. Local `.h5ad` fallback (any `.h5ad` in `/tmp` or `work/`).
4. Pure-R synthetic fallback (default; no external dependencies).

The synthetic fallback builds a small Seurat-like list object with:
- 180 cells, 500 HVGs (RNA), 2000 peaks (ATAC).
- 3 cell types (`CT1`, `CT2`, `CT3`) for the `CellType` column.
- 10 TF symbols (GATA1, MYB, KLF1, ...) injected into the RNA assay.
- UMAP.rna and UMAP.atac reductions, 2-D each.

## Use with the Nextflow pipeline

```bash
# Smoke test on a laptop without the full 6.5 GB raw data:
nextflow run main.nf \
    --input  data_demo/synthetic/mini_seurat.rds \
    --percentile 90 \
    --execution_mode synthetic \
    -profile test
```

## Why a synthetic object at all?

CellOracle and SeuratDisk have a slow first-run cold start (R package
compilation, motif downloads). The synthetic object lets the CI verify the
**structural** correctness of the pipeline (channel binding, output
schemas, file publication) in under 30 s without paying the cost of the
biological dataset.

The biological object produced by `data_demo/subsample_dataset.py` is the
one that produces real, scientifically meaningful CellOracle output. Use the
synthetic object only for pipeline plumbing tests.
