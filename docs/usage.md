# Usage

## Prerequisites

- [Nextflow](https://www.nextflow.io/) >= 25.10.0
- [Docker](https://www.docker.com/), [Singularity](https://docs.sylabs.io/guides/latest/user_guide/) or
  [Conda](https://docs.conda.io/) (one of them)
- Java >= 11

## Run the pipeline

```bash
# 1) (one-time) obtain and subsample the biological test dataset
bash  data_demo/download_data.sh        # ~30 min, ~6.5 GB into data_demo/raw/
python data_demo/subsample_dataset.py   # ~10 min, produces data_demo/processed/

# 2) Run with the biological microglia dataset and Docker
nextflow run main.nf -profile test,docker

# 3) Re-run with a different network-pruning percentile.
#    The heavy R preprocessing step is CACHED; only Python re-runs.
nextflow run main.nf -profile test,docker --percentile 95 -resume
```

## Profiles

| Profile   | Purpose                                                         |
|-----------|-----------------------------------------------------------------|
| `test`    | Use the bundled `data_demo/processed/microglia_dam_demo.rds` (the biological microglia dataset). |
| `docker`  | Run every task in a Docker container.                           |
| `singularity` | Run every task in a Singularity container (HPC-friendly).   |
| `conda`   | Use the unified `containers/environment.yml`.                   |
| `standard` | Default: rely on whatever the local environment provides.     |

## Key parameters

| Parameter      | Default                                                       | Description                                              |
|----------------|---------------------------------------------------------------|----------------------------------------------------------|
| `--input`      | `data_demo/processed/microglia_dam_demo.rds`                  | Input Seurat RDS object.                                 |
| `--outdir`     | `results`                                | Output directory.                                        |
| `--percentile` | `90`                                     | Network link pruning percentile (the `-resume` demo param). |
| `--genome`     | `hg38`                                   | Reference genome assembly.                               |
| `--ko_targets` | `GATA1,MYB,KLF1`                         | Comma-separated TF symbols for in silico KO simulation.  |
