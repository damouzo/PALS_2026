# Troubleshooting

## Common issues

### `command not found: nextflow`

Install Nextflow:

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/   # or any directory on $PATH
```

### `OutOfMemoryError` in `INFER_GRN`

CellOracle's `knn_imputation` step is memory-intensive. Increase the
`process_high` label in `conf/base.config`:

```groovy
withLabel: 'process_high' {
    memory = '64.GB'
}
```

### `Docker` daemon not running

```bash
# macOS
open -a Docker
# Linux
sudo systemctl start docker
```

### `celloracle` cannot download motif database

The pipeline expects an offline motif database. If `motif_db` cannot be fetched
at runtime, pre-cache it in the Docker image (Phase 7) and reference its path
via the `--motif_db` parameter.

### Re-running after editing the R script

If you change `bin/GRN_dataProcess.R`, Nextflow will detect a different task
hash and re-run `PREPROCESS_ATAC` even with `-resume`. This is by design.

### `nextflow clean` to wipe the cache

```bash
# Remove all task work directories (force)
nextflow clean -f
```
