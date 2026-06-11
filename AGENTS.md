This file contains everything an AI agent (or new developer) needs to understand and work on this codebase effectively.

> **Communication Language**: The user prefers to communicate in **Spanish**. Always respond in Spanish when interacting with the user, unless they explicitly switch languages. However, all code, documentation, logs, and commit messages in this repository must remain in **English**.

---

## Project Summary

**scMultiome-GRN** is a Nextflow DSL2 pipeline designed to integrate Single-Cell Multiome data (scRNA-seq + scATAC-seq), infer Gene Regulatory Networks (GRNs), and simulate *in silico* cellular perturbations. 

The project features a native polyglot design: it leverages **R (Seurat/Signac/Cicero)** for the initial chromatin accessibility preprocessing and co-accessibility analysis, and **Python (Scanpy/CellOracle)** for predictive network inference and dynamic cellular state simulations.

It is surgically optimized for a **1-hour Live-Demo / Masterclass**. It runs an ultra-downsampled public dataset on the presenter's laptop in under 2 minutes, spectacularly demonstrating the power of Nextflow's caching system (`-resume`).

---

## Repository Layout

scMultiome-GRN/
├── bin/                     # Clean executable scripts called by Nextflow processes
│   ├── GRN_dataProcess.R    # ATAC/RNA preprocessing & Cicero co-accessibility integration
│   └── GRN_analysis.py      # Network inference with CellOracle & perturbation simulations
├── conf/                    # Nextflow configuration files
│   ├── base.config          # Default resource allocation profiles
│   └── modules.config       # Process-specific publishDir directives and CLI argument mappings
├── containers/
│   ├── Dockerfile           # PRIMARY — Unified R + Python environment (used by GitHub Actions/Docker)
│   └── environment.yml      # Unified Conda environment definition for local development/fallback
├── data_demo/               # Everything related to test data: download, subsample, docs
│   ├── raw/                 # Raw 10x assets (gitignored, ~6.5 GB)
│   ├── processed/           # Subsampled outputs (gitignored, generated)
│   ├── download_data.sh     # Idempotent curl-based downloader
│   ├── subsample_dataset.py # Scanpy-based subsampler (microglia from 10x clusters)
│   ├── README.md            # Biological + technical documentation
│   └── synthetic/           # Optional ultra-tiny synthetic object (smoke tests)
├── modules/local/           # Nextflow DSL2 process definitions (one directory per process)
│   ├── preprocess_atac.nf   # Wraps the R preprocessing script execution
│   └── infer_grn.nf         # Wraps the Python network inference script execution
├── workflows/
│   └── multiome_grn.nf      # Workflow orchestration and channel binding (R -> Python)
├── main.nf                  # Pipeline entry point
└── nextflow.config          # Primary Nextflow configuration (imports conf/)


---

## Key Design Decisions

### 1. Elimination of Legacy Hacks (No YAML Parsing in Bash)
The old `GRN.sh` script utilized an external shell script (`parse_yaml.sh`) to inject static variables from a YAML file into Bash environments. This design has been completely deprecated. Variables are now natively managed via Nextflow `params` and passed to the scripts in `bin/` as clean, structured command-line arguments using standard flags (`--input`, `--outdir`, `--percentile`).

### 2. Unified Polyglot Environment (Single Container)
Because the pipeline transitions sequentially from R to Python, `containers/Dockerfile` and `containers/environment.yml` act as the single source of truth. They package both R genomic libraries (`Seurat`, `Signac`, `cicero`) and the Python ecosystem (`scanpy`, `celloracle`) together. This ensures seamless cross-language task execution without environment or container switching overhead.

### 3. Parameter Isolation for Caching Demonstration (`-resume`)
The `--percentile` parameter (used by CellOracle to prune weak network links) is declared as an isolated top-level parameter in Nextflow (`params.percentile`). This enables the presenter to run the full pipeline once, modify this single numerical value on the fly, and re-run using the `-resume` flag. Nextflow detects that the heavy R preprocessing step (`PREPROCESS_ATAC`) remains untouched, uses the cache instantaneously, and executes only the Python step (`INFER_GRN`) in seconds in front of the audience.

### 4. Dictatorial Test Dataset Scale
The biological test dataset (`data_demo/processed/microglia_dam_demo.rds`) must be strictly constrained in dimensions (~500–1000 cells, restricted to a few chromosomes, e.g. chr7 + chr15 + chr17). Any future updates must validate that Cicero takes less than 1 minute to compute co-accessibilities on standard dual-core laptop hardware. The optional `data_demo/synthetic/` ultra-tiny object (~ 300 KB, 150–200 cells) is reserved for CI smoke tests where the full biological dataset is unavailable.

---

## Container & CI/CD

### Source of Truth
The **Dockerfile** (`containers/Dockerfile`) is the authoritative environment definition:
- Integrates all required R bioinformatic packages and Python analytical frameworks.
- Automated via **GitHub Actions**: The image builds and pushes automatically to GitHub Container Registry (GHCR) at `ghcr.io/dmouzo/scmultiome-grn` on every push to `main`, on `pull_request` (build only, no push), and on `v*.*.*` tags.
- The `environment.yml` file is maintained purely for local development flexibility to bootstrap the environment on a laptop without launching Docker.

### Pinning Strategy (validated June 2026)
- `r-base=4.4.3` (r44 build-string; the only base version that resolves
  against the modern r-seurat / r-signac stack on conda-forge).
- `r-seurat=5.1.0` INCLUDES SeuratWrappers since Seurat v5. The legacy
  `r-seuratwrappers` package is abandoned on conda-forge; do NOT add it
  back. This was the original cause of the broken build in Phase 7.
- `r-signac=1.14.0` + `r-ggnewscale` (the latter is required for
  modern coverage plots and was missing in the legacy image).
- `r-seuratdisk` is NOT installable from conda-forge for the r44
  build-string (0.0.9019 dates from 2021 and depends on r-spatstat<2.0
  which has no r44-compatible build). It is installed from CRAN in
  the Dockerfile post-conda RUN. `GRN_dataProcess.R` has a fallback
  to .rds sidecars if SeuratDisk is missing.
- `r-monocle3` and `r-cicero` are Bioconductor-only; installed via
  `BiocManager::install(..., version = "3.19")` post-conda.
- `celloracle==0.18.0` + `gimmemotifs==3.0.0` mirror the upstream
  `kenjikamimoto126/celloracle_ubuntu:0.18.0` image (used as a pin
  reference, not as a base).

### Local Development
- Building the full image locally requires ≥8 GB RAM and ~20 min on a
  modern laptop. The `mambaorg/micromamba:1.5.10-jammy` base image is
  used for fast R+Python installation via `micromamba install`.
- The Dockerfile ends with a `RUN` that imports every required R and
  Python package. If any dependency is missing the image WILL NOT be
  pushed to GHCR — the GitHub Actions job fails before the push step.
- For local-only iteration (no Docker), bootstrap a conda env from
  `environment.yml` with `mamba env create -f containers/environment.yml`
  then activate it and run Nextflow with `-profile conda`.

---

## Running the Pipeline (Live Demo Commands)

```bash
# PRE-DEMO: obtain and prepare the biological test dataset
bash data_demo/download_data.sh        # one-time, ~30 min, ~6.5 GB
python data_demo/subsample_dataset.py  # ~10 min, produces data_demo/processed/

# LIVE DEMO - STEP 1: Full initial execution (R preprocessing + Python inference)
nextflow run main.nf -profile test,docker

# LIVE DEMO - STEP 2: The Cache WOW-Factor
# Modify the network pruning percentile on the fly.
# The heavy R step is skipped completely (CACHED), running only Python in seconds.
nextflow run main.nf -profile test,docker --percentile 95 -resume

Development Conventions
Nextflow Workflow & DSL2

    Every process must be encapsulated inside its own module under modules/local/.

    Channels must connect the output of R (all_peaks.csv and cicero_connections.csv) cleanly into the input channels of the Python process.

    Nextflow processes must be named in UPPER_SNAKE_CASE.

    Nextflow channels must be named in ch_camelCase.

R Scripts (bin/GRN_dataProcess.R)

    Must use optparse or argparse to cleanly capture --input_seurat and --outdir.

    Must clear Seurat objects from memory after exporting structural data to prevent Out-Of-Memory (OOM) errors on limited local laptop hardware.

    Must generate deterministic, unquoted CSV tables required as inputs for the Python step.

Python Scripts (bin/GRN_analysis.py)

    Must use argparse to capture individual paths for CSV inputs and the numerical --percentile value.

    All final plots generated by celloracle (virtual knock-out vector fields and cell fate shift plots) must be saved dynamically using the path provided by the --outdir parameter.

Workflow Data Flow

     mini_seurat.rds (Input Channel)
               ↓
     [ PREPROCESS_ATAC ]  <--- (R Process: Signac / Cicero)
        /           \
 all_peaks.csv     cicero_connections.csv  (Data Channels)
        \           /
       [ INFER_GRN ]      <--- (Python Process: CellOracle / params.percentile)
               ↓
     Output Plots & Figures (PublishDir)

Known Limitations

    Optimized for Local Demo Only: The pipeline defaults assume highly restricted memory and CPU constraints suitable for quick local execution. Scaling this to large patient multiomic cohorts on HPC clusters will require allocating specific high-resource limits inside conf/base.config.

    Internet Dependency for Motifs: Certain internal routines within celloracle download Transcription Factor motif databases on demand. The developer agent must ensure that standard base motif dictionaries (or downsampled reference motifs) are pre-cached or packaged within the environment to avoid relying on the event venue's live Wi-Fi connection.