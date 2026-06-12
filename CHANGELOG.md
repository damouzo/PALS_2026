# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Multi-source container architecture** replacing the previous
  GHCR + Wave mix:
  - `pals-r-seurat` (Wave ORAS link
    `oras://community.wave.seqera.io/library/python_r-base_r-essentials_r-devtools_pruned:79376aa02ef4e2df`):
    R 4.4.3 + R-essentials + R-devtools + tidyverse + r-seurat 5.4.0 +
    hdf5r + optparse + the Bioconductor base classes. Signac and
    SeuratDisk are installed at runtime by `bin/GRN_seurat.R`. Used by
    `PREP_R_SEURAT`.
  - `pals-r-bioc` (Wave ORAS link
    `oras://community.wave.seqera.io/library/python_r-base_r-essentials_r-remotes_pruned:ffe857cbb37c1f02`):
    R 4.4.3 + R-essentials + R-remotes + the Bioconductor base classes.
    cicero + monocle3 + EnsDb.Mmusculus.v79 are installed at runtime by
    `bin/GRN_cicero.R`. Used by `PREP_R_CICERO`.
  - `pals-python-grn:1.0.0` (`cr.seqera.io/dmouzo/`): re-tag of
    `kenjikamimoto126/celloracle_ubuntu:0.18.0` (CellOracle 0.18.0 +
    GimmeMotifs 0.18.3 + scanpy + anndata + goatools + adjustText +
    scientific stack). Used by `INFER_GRN`.
- **Sequera Wave sources** under `containers/`:
  - `containers/env-r-seurat.yml` — historical conda spec for the seurat
    image (marked LEGACY; the Wave ORAS link is the live source).
  - `containers/env-r-bioc.yml` — historical conda spec for the cicero
    image (marked LEGACY).
  - `containers/Dockerfile.r-seurat` — historical post-conda Dockerfile
    (marked LEGACY).
  - `containers/Dockerfile.r-bioc` — historical post-conda Dockerfile
    (marked LEGACY).
- **New Nextflow parameters**: `container_r1`, `container_r2`,
  `container_py` (overridable per run). Legacy `container_r` alias is
  kept and maps to `container_r1` for back-compat.
- `data_demo/` ecosystem for reproducible biological test data:
  - `data_demo/raw/` — 10x Genomics raw assets (gitignored, ~6.5 GB).
  - `data_demo/processed/` — Subsampled outputs (gitignored, generated).
  - `data_demo/download_data.sh` — Idempotent curl-based downloader (4 assets, ~30 min).
  - `data_demo/subsample_dataset.py` — Scanpy-based subsampler that uses the OFFICIAL
    10x Cell Ranger ARC 2.0.1 cluster labels (GEX cluster 21 = microglia,
    cluster 22 = DAM-like) to extract 3 biologically meaningful CellTypes:
    `Microglia_Homeostatic` (WT 13.2m), `Microglia_EarlyDAM` (TgCRND8 5.7m),
    `Microglia_LateDAM` (TgCRND8 17.9m). Symmetric 200 cells per CellType.
    - ATAC pruned to chr7 + chr15 + chr17 (Apoe / Slc1a3 / Trem2 loci).
    - Co-accessibility input includes Slc1a3 distal peak highlighted in the 10x AppNote.
    - Tabix-accelerated fragment filtering (with pure-Python fallback).
    - Lazy imports: `--help` works without scanpy/anndata installed.
    - CLI supports a fully customisable `--cohorts` YAML/JSON spec.
  - `data_demo/synthetic/` — Ultra-tiny synthetic Seurat object (~300 KB) reserved
    for CI smoke tests where the multi-GB raw data is unavailable.
  - `data_demo/.gitignore` (defence in depth).
  - `data_demo/README.md` — Biological + technical documentation, references.

### Changed
- **Container source of truth is now a Wave + GHCR mix**: the two R-side
  images are Sequera Wave **ORAS links** (the community cache resolves
  them on the executor) and the Python-side image is re-tagged to
  **GHCR** under `ghcr.io/damouzo/pals-python-grn:1.0.0` (note the
  GitHub handle is `damouzo`, not `dmouzo`). The earlier
  `ghcr.io/dmouzo/...` and `cr.seqera.io/dmouzo/pals-r-*` references
  have been removed from `nextflow.config`,
  `nextflow_schema.json`, `AGENTS.md`, `README.md` and
  `data_demo/subsample_dataset.py`. Signac, SeuratDisk, cicero,
  monocle3 and EnsDb.Mmusculus.v79 are no longer baked into the
  R-side images; they are installed at runtime by `bin/GRN_seurat.R`
  and `bin/GRN_cicero.R` via idempotent helpers.
- **CI removed**: `.github/workflows/ci-tests.yml` deleted. Validation
  is now manual via `nextflow lint .` plus a local
  `nextflow run main.nf -profile test,docker` smoke test.
- **Pipeline architecture**: R-side preprocessing is now split into TWO
  specialized processes (`PREP_R_SEURAT` then `PREP_R_CICERO`) instead of
  a single monolithic `PREPROCESS_ATAC`. This allows each process to run
  in its own dedicated image so the conda solver stays tractable.
  - `bin/GRN_dataProcess.R` → split into `bin/GRN_seurat.R` +
    `bin/GRN_cicero.R`.
  - `modules/local/preprocess_atac.nf` → replaced by
    `modules/local/prep_r_seurat.nf` + `modules/local/prep_r_cicero.nf`.
  - `workflows/multiome_grn.nf` updated to chain the three processes.
- **Architecture**: `test_data/` removed. All test-data assets now live in
  `data_demo/`. The default pipeline `--input` points to
  `data_demo/processed/microglia_dam_demo.rds`.
- **Environment**: `containers/environment.yml` switched
  `bioconductor-ensdb.hsapiens.v86` → `bioconductor-ensdb.mmusculus.v79`
  (Ensembl 79 / GRCm38) to match the 10x Mouse Brain dataset.
- `bin/GRN_cicero.R` now uses `EnsDb.Mmusculus.v79` (the previous
  hardcoded `EnsDb.Hsapiens.v86` was a bug introduced during the
  dataset migration and is fixed in this release).
- `nextflow.config`, `nextflow_schema.json`, `AGENTS.md`, `README.md`,
  `docs/parameters.md`, `docs/usage.md`, `.github/workflows/ci-tests.yml`
  all updated to point to the new biological default input and the new
  three-container architecture.

### Fixed
- **Container architecture**: replaced the self-hosted Docker image
  (`containers/Dockerfile` + `containers/environment.yml` built by
  `.github/workflows/docker-publish.yml` and pushed to GHCR) with **two
  pre-built images hosted on Seqera Containers**, one per process:
  - `cr.seqera.io/dmouzo/pals-r-preprocess:1.0.0` (R 4.4.3 + Seurat
    5.1.0 + Signac 1.14.0 + SeuratDisk + cicero + monocle3 +
    EnsDb.Mmusculus.v79 + utility R packages).
  - `cr.seqera.io/dmouzo/pals-python-grn:1.0.0` (Python 3.11 +
    celloracle 0.18.0 + scanpy + anndata + gimmemotifs 0.18.3 +
    goatools + adjustText + scientific stack).
- **New parameters**: `params.container_r` and `params.container_py`
  (overridable per run with `--container_r <url> --container_py <url>`).
  Both default to the Seqera Containers URLs.
- **Removed legacy GHA workflow**: `.github/workflows/docker-publish.yml`
  is gone (we no longer build images in CI).
- **Simplified `.github/workflows/ci-tests.yml`**: the smoke test now
  pulls the two pre-built Seqera Containers via `docker run` for the
  R-only data prep / validation steps, and the pipeline run uses them
  via `-profile docker`. Build times drop from ~30 min to ~3 min on CI.
- **`gimmemotifs` packaging**: the source distribution on PyPI is broken
  on Python 3.11+ (depends on `configparser.SafeConfigParser` removed in
  3.12). We pin `gimmemotifs=0.18.3` from bioconda, where the upstream
  maintainer publishes pre-built py311 wheels.
    SeuratDisk is now installed from CRAN in the Dockerfile post-conda
    RUN, with the GRN_dataProcess.R fallback to .rds sidecars if it
    cannot be installed.
  - `r-monocle3` and `r-cicero` (Bioconductor-only) installed via
    `BiocManager::install(..., version = "3.19")` post-conda.
  - `celloracle==0.18.0` + `gimmemotifs==3.0.0` aligned with the upstream
    `kenjikamimoto126/celloracle_ubuntu:0.18.0` image (used as a pin
    reference, not as a base).
  - Added `pandoc`, `libfontconfig1`, `libxt-dev` to the apt layer.
- **Dockerfile**: added a post-build validation `RUN` that imports every R
  and Python package. If any dependency is missing the image WILL NOT
  be pushed to GHCR. Pinned `LANG=C.UTF-8`, `PYTHONNOUSERSITE=1`,
  `RETICULATE_PYTHON=/opt/conda/bin/python3`, `MAMBA_NO_LOCKFILE_INSTALL=1`
  for deterministic runtime.
- **CI workflows**:
  - `docker-publish.yml`: multi-arch (linux/amd64 + linux/arm64), explicit
    post-build R + Python verification, deterministic tag (sha-${{ github.sha_short }}).
  - `ci-tests.yml`: smoke test now runs INSIDE the published GHCR image
    (needs: docker-publish) instead of installing micromamba on the runner.
    Wall-clock time for the smoke job drops from ~6 min to ~30 s on
    subsequent runs (GHCR pull is fast).

### Removed
- `containers/` directory (Dockerfile + environment.yml) and
  `.github/workflows/docker-publish.yml`. Container images are no longer
  built in CI; they are pulled from Seqera Containers at runtime.
- `params.container_tag` (replaced by `params.container_r` and
  `params.container_py`, both with full image URLs as defaults).

## [0.1.0] - 2024-06-09

### Added
- Initial Nextflow DSL2 pipeline skeleton (`main.nf`, `nextflow.config`).
- DSL2 process modules: `PREPROCESS_ATAC` (R/Seurat/Signac/Cicero) and `INFER_GRN` (Python/CellOracle).
- Reusable sub-workflow (`workflows/multiome_grn.nf`) that orchestrates the R -> Python handoff.
- `bin/GRN_dataProcess.R` refactored as a clean CLI tool (`optparse` + manual `commandArgs` fallback for CI).
  - Supports full Seurat + Signac + Cicero pipeline when those packages are available.
  - Includes a deterministic synthetic-fallback path for CI smoke tests.
  - Exports `all_peaks.csv`, `cicero_connections.csv`, `rna.h5ad`, `atac.h5ad`.
  - Memory-hygiene pattern: explicit `rm(); gc()` before exit.
- `bin/GRN_analysis.py` refactored as a clean CLI tool (`argparse`).
  - Two execution modes: `full` (CellOracle + GimmeMotifs) and `synthetic` (offline CI fallback).
  - Auto-detection via `mode=auto` (default).
  - In silico perturbation module with vector-field plots for knockout targets (GATA1, MYB, KLF1 by default).
  - Resilient serialization: parquet when `pyarrow`/`fastparquet` is available, CSV fallback otherwise.
  - Excel/CSV dual output for `network_scores`.
- `bin/build_test_data.R` for regenerating the synthetic test dataset.
- `bin/validate_test_data.R` for sanity-checking the test dataset.
- Unified polyglot `Dockerfile` (micromamba + R 4.3 + Python 3.11 + Java 17 + Nextflow).
- Fallback `environment.yml` (conda-forge + bioconda).
- GitHub Actions workflows:
  - `docker-publish.yml`: builds and pushes the image to GHCR on push to main, PRs, and releases.
  - `ci-tests.yml`: lint + smoke test on every CI run.
- Test dataset (`test_data/mini_seurat.rds`): 180 cells, 500 HVGs, 2000 peaks, 3 cell types, 0.29 MB.
- Five profiles: `standard`, `test`, `docker`, `singularity`, `conda`.
- Per-label resource allocation (`process_low`, `process_medium`, `process_high`) overridable per profile.
- `nextflow_schema.json` with 15 documented parameters + enums + defaults.
- Automatic reports: `execution_timeline.html`, `execution_report.html`, `execution_trace.txt`, `pipeline_dag.svg`.
- Documentation: `docs/usage.md`, `docs/output.md`, `docs/parameters.md`, `docs/troubleshooting.md`.

### Changed
- Refactored `docs/GRN/GRN.sh` legacy Bash wrapper out of the active pipeline.
- Replaced hardcoded paths and external `parse_yaml.sh` dependency with `params`.
- Eliminated duplicated code blocks in the legacy `GRN_analysis.py`.
- Dynamic tag in `INFER_GRN` (`p${percentile}_${ko_targets}`) makes the `-resume` cache demo immediately visible in the log.

### Removed
- `docs/GRN/GRN.sh` (replaced by Nextflow orchestration; kept only as historical reference).
- External `parse_yaml.sh` Bash dependency.

### Fixed
- Output glob patterns in `INFER_GRN` accept both `.parquet` and `.csv` (or `.xlsx`/`.csv`) so the pipeline works regardless of which serialization libs are present in the runtime.

### Security
- Pipeline runs as the non-root `mamaba_user` inside the container.
- Minimal CI permissions: `packages: write`, `contents: read`.
