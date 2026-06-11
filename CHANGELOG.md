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
- **Architecture**: `test_data/` removed. All test-data assets now live in
  `data_demo/`. The default pipeline `--input` points to
  `data_demo/processed/microglia_dam_demo.rds`.
- **Environment**: `containers/environment.yml` switched
  `bioconductor-ensdb.hsapiens.v86` → `bioconductor-ensdb.mmusculus.v79`
  (Ensembl 79 / GRCm38) to match the 10x Mouse Brain dataset.
- `bin/GRN_dataProcess.R` already uses `EnsDb.Mmusculus.v79` if the container
  is built from the updated `environment.yml`.
- `nextflow.config`, `nextflow_schema.json`, `AGENTS.md`, `README.md`,
  `docs/parameters.md`, `docs/usage.md`, `.github/workflows/ci-tests.yml`
  all updated to point to the new biological default input.

### Fixed
- **Container build (`containers/environment.yml` + `Dockerfile`)**:
  rewritten from scratch to fix the long-standing resolver failure. Pins
  now validated against conda-forge/bioconda (June 2026):
  - `r-base=4.4.3` (r44 build-string; only version that resolves against
    the modern r-seurat / r-signac stack).
  - `r-seurat=5.1.0` (Seurat v5 INCLUDES SeuratWrappers; the legacy
    `r-seuratwrappers` package is abandoned on conda-forge).
  - `r-signac=1.14.0` + `r-ggnewscale` (required for modern coverage plots).
  - `r-seuratdisk=0.0.9019` REMOVED from conda deps (broken against r44;
    depends on `r-spatstat<2.0` which has no r44-compatible build).
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
- `containers/test_minimal.yml` (unreferenced residue from earlier
  iteration; not used by any workflow or script).

### Removed
- `test_data/` (the synthetic 180-cell mini object) and the scripts that
  produced/validated it (`bin/build_test_data.R`, `bin/validate_test_data.R`).
  Their functionality lives on in `data_demo/synthetic/build_synthetic_mini_seurat.R`.
- `docs/GRN/` directory containing the legacy `GRN.sh`, `GRN_analysis.py`,
  `GRN_dataProcess.R` and WSL `Zone.Identifier` artifacts.
- `containers/test_minimal.yml` (see [Unreleased] > Fixed above).

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
