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
│   ├── GRN_seurat.R         # Stage 1 R: scRNA + scATAC preprocessing (Seurat + Signac + SeuratDisk)
│   ├── GRN_cicero.R         # Stage 2 R: Cicero co-accessibility (cicero + monocle3 + EnsDb.Mmusculus.v79)
│   └── GRN_analysis.py      # Stage 3 Python: Network inference with CellOracle & perturbation simulations
├── conf/                    # Nextflow configuration files
│   ├── base.config          # Default resource allocation profiles
│   └── modules.config       # Process-specific publishDir directives and CLI argument mappings
├── containers/              # Sequera Wave sources for the R-side images
│   ├── env-r-seurat.yml     # Conda spec for pals-r-seurat (R 4.4.3 + Seurat + Signac + SeuratDisk)
│   ├── env-r-bioc.yml       # Conda spec for pals-r-bioc (R 4.4.3 + BiocManager → cicero/monocle3/EnsDb)
│   ├── Dockerfile.r-seurat  # Post-conda installs (BiocManager, Signac, SeuratDisk) for pals-r-seurat
│   └── Dockerfile.r-bioc    # Post-conda installs (BiocManager → cicero/monocle3/EnsDb) for pals-r-bioc
├── data_demo/               # Everything related to test data: download, subsample, docs
│   ├── raw/                 # Raw 10x assets (gitignored, ~6.5 GB)
│   ├── processed/           # Subsampled outputs (gitignored, generated)
│   ├── download_data.sh     # Idempotent curl-based downloader
│   ├── subsample_dataset.py # Scanpy-based subsampler (microglia from 10x clusters)
│   ├── README.md            # Biological + technical documentation
│   └── synthetic/           # Optional ultra-tiny synthetic object (smoke tests)
├── modules/local/           # Nextflow DSL2 process definitions (one file per process)
│   ├── prep_r_seurat.nf     # PREP_R_SEURAT: R 4.4.3 + Seurat + Signac + SeuratDisk
│   ├── prep_r_cicero.nf     # PREP_R_CICERO: R 4.4.3 + cicero + monocle3 + EnsDb.Mmusculus.v79
│   └── infer_grn.nf         # INFER_GRN: Python 3.11 + CellOracle + GimmeMotifs
├── workflows/
│   └── multiome_grn.nf      # Workflow orchestration: R -> R -> Python chain
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
The pipeline uses **three pre-built container images**, one per process.
There is no in-repo Dockerfile and no in-CI image build step. The R-side
images are built with **Sequera Wave** (UI-driven conda solver on beefy
infrastructure); the Python-side image is a re-tag of
`kenjikamimoto126/celloracle_ubuntu:0.18.0` published on **GHCR**.

- **`pals-r-seurat`** (`cr.seqera.io/dmouzo/pals-r-seurat:1.0.0`):
  R 4.4.3 + Seurat 5.4.0 + Signac 1.14.0 + SeuratDisk 0.0.9019 +
  utility R packages (tidyverse, matrix, rccp*, hdf5r, optparse) +
  Bioconductor base (Biobase, S4Vectors, IRanges, GenomicRanges, etc.).
  Used by `PREP_R_SEURAT`.
- **`pals-r-bioc`** (`cr.seqera.io/dmouzo/pals-r-bioc:1.0.0`):
  R 4.4.3 + minimal Bioconductor base + `BiocManager::install` →
  cicero + monocle3 + EnsDb.Mmusculus.v79. Used by `PREP_R_CICERO`.
- **`pals-python-grn`** (`ghcr.io/dmouzo/pals-python-grn:1.0.0`):
  Python 3.11 + celloracle 0.18.0 + scanpy + anndata +
  gimmemotifs 0.18.3 + goatools + adjustText + scientific stack.
  Used by `INFER_GRN`.

The three image URLs are exposed as top-level Nextflow params
(`params.container_r1`, `params.container_r2`, `params.container_py`)
and can be overridden per run with `--container_r1 <url>`,
`--container_r2 <url>`, `--container_py <url>`. The legacy alias
`--container_r` is mapped to `--container_r1` for back-compat.
The `nextflow_schema.json` documents every parameter.

### Why three images (and not one)?
After three failed CI builds (each surfacing a different solver or
pinning issue: `r-seuratdisk` requires an older `r-spatstat`,
`gimmemotifs` on PyPI is broken on Python 3.11+, the joint
`r-seurat + r-signac + monocle3 + EnsDb` solver does not converge on
libmamba), we split the R-side environment into two specialized images
so the conda solver stays tractable. The Python-side environment is a
re-tag of the upstream `kenjikamimoto126/celloracle_ubuntu:0.18.0` image
so we never re-derive the (extremely fragile) `gimmemotifs`-on-Py3.11
installation ourselves.

### Why Sequera Containers (and not a self-built image in CI)?
The conda solver for the R stack requires a beefier environment than
the GitHub Actions runner provides. Sequera Wave runs the solver on
its own infrastructure and the build is inspectable step-by-step. The
CI workflow now just `docker pull`s the three images and runs the
pipeline.

### CI
- `.github/workflows/ci-tests.yml`:
  - `lint` job: `nextflow lint .` on every push.
  - `smoke` job: pulls the three containers (2 from Sequera Wave, 1
    from GHCR), builds the synthetic test dataset inside the
    R-seurat image, then runs the full pipeline via
    `nextflow run main.nf -profile test,docker`.
- There is no `docker-publish.yml` workflow (we no longer build images
  in CI).

### Local Development
- For local execution with the three images, no setup is needed beyond
  `docker pull` and a working `nextflow` install:
  ```bash
  docker pull cr.seqera.io/dmouzo/pals-r-seurat:1.0.0
  docker pull cr.seqera.io/dmouzo/pals-r-bioc:1.0.0
  docker pull ghcr.io/dmouzo/pals-python-grn:1.0.0
  nextflow run main.nf -profile test,docker
  ```
- For local iteration on a single process, override the container URL
  with a local tag:
  ```bash
  nextflow run main.nf -profile test,docker \
      --container_r1 'docker://my-r-seurat:dev' \
      --container_r2 'docker://my-r-bioc:dev' \
      --container_py 'docker://my-python-grn:dev'
  ```
- The `standard` profile (no Docker) requires R 4.4.x and Python 3.11
  installed on the host with the same packages. Not recommended for
  the demo — the Sequera + GHCR images are the supported runtime.

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

R Scripts (bin/GRN_seurat.R, bin/GRN_cicero.R)

    Must use optparse or argparse to cleanly capture their CLI flags.

    Must clear Seurat objects from memory after exporting structural data to prevent Out-Of-Memory (OOM) errors on limited local laptop hardware.

    Must generate deterministic, unquoted CSV tables required as inputs for the Python step.

    GRN_cicero.R must use `EnsDb.Mmusculus.v79` (NOT `EnsDb.Hsapiens.v86`) since the demo dataset is mm10/mouse.

Python Scripts (bin/GRN_analysis.py)

    Must use argparse to capture individual paths for CSV inputs and the numerical --percentile value.

    All final plots generated by celloracle (virtual knock-out vector fields and cell fate shift plots) must be saved dynamically using the path provided by the --outdir parameter.

Workflow Data Flow

     mini_seurat.rds (Input Channel)
               ↓
     [ PREP_R_SEURAT ]   <--- (R Process: Seurat + Signac + SeuratDisk)
        /        |        \
   rna.h5ad  atac.h5ad  seurat_object.rds
        \        |        /
     [ PREP_R_CICERO ]   <--- (R Process: cicero + monocle3 + EnsDb.Mmusculus.v79)
        /           \
 all_peaks.csv     cicero_connections.csv
        \           /
       [ INFER_GRN ]      <--- (Python Process: CellOracle / params.percentile)
               ↓
     Output Plots & Figures (PublishDir)

Known Limitations

    Optimized for Local Demo Only: The pipeline defaults assume highly restricted memory and CPU constraints suitable for quick local execution. Scaling this to large patient multiomic cohorts on HPC clusters will require allocating specific high-resource limits inside conf/base.config.

    Internet Dependency for Motifs: Certain internal routines within celloracle download Transcription Factor motif databases on demand. The developer agent must ensure that standard base motif dictionaries (or downsampled reference motifs) are pre-cached or packaged within the environment to avoid relying on the event venue's live Wi-Fi connection.