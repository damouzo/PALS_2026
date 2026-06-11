# Parameters

All parameters are declared in `nextflow.config` and documented in
`nextflow_schema.json`. The table below summarizes the most important ones.

| Parameter      | Default                                  | Description                                                                          |
|----------------|------------------------------------------|--------------------------------------------------------------------------------------|
| `--input`      | `data_demo/processed/microglia_dam_demo.rds` | Path to the input Seurat RDS object (the biological microglia dataset).            |
| `--outdir`     | `results`                                | Output directory.                                                                    |
| `--genome`     | `hg38`                                   | Reference genome assembly.                                                           |
| `--motif_db`   | `gimme.vertebrate.v5.0`                  | GimmeMotifs database identifier.                                                     |
| `--clustering` | `CellType`                               | `meta.data` column holding the cell-type clustering.                                 |
| `--umap_key`   | `umap_rna`                               | Reduction key in the Seurat object.                                                  |
| `--clust2study`| `HSC_1`                                  | Cluster used as the anchor for in silico perturbation.                               |
| `--percentile` | `90`                                     | Network link pruning percentile. **Isolated** as a top-level param to showcase `-resume`. |
| `--p_threshold`| `0.001`                                  | P-value threshold for link filtering.                                                |
| `--n_top_edges`| `2000`                                   | Maximum number of edges kept after filtering.                                        |
| `--ko_targets` | `GATA1,MYB,KLF1`                         | Comma-separated TF symbols to simulate as knockouts.                                 |
| `--max_cpus`   | `8`                                      | Maximum CPUs that any single task can request.                                       |
| `--max_memory` | `32.GB`                                  | Maximum memory that any single task can request.                                     |
| `--max_time`   | `4.h`                                    | Maximum wall-clock time for any single task.                                         |
| `--container_tag` | `latest`                              | Docker/Singularity image tag.                                                        |

## The `--percentile` parameter: the `-resume` demo driver

`--percentile` is the only parameter that controls the **Python** post-processing
step. The R preprocessing step is **insensitive** to it. This is the key
architectural choice that makes the live-demo `-resume` trick work:

```bash
# First run: full execution (R + Python)
nextflow run main.nf -profile test,docker

# Change only the percentile; -resume re-uses the R cache and re-runs Python only.
nextflow run main.nf -profile test,docker --percentile 95 -resume
```

> The default `--input` is the biological microglia dataset. For CI smoke tests
> where the multi-GB raw data is unavailable, point `--input` to
> `data_demo/synthetic/mini_seurat.rds` (an ultra-tiny synthetic object regenerated
> by `data_demo/synthetic/build_synthetic_mini_seurat.R`).
