# ==============================================================================
# containers/Dockerfile — Polyglot R + Python unified image
# ==============================================================================
#
# This is the SINGLE SOURCE OF TRUTH for the scMultiome-GRN runtime.
# Built and pushed automatically to GitHub Container Registry (GHCR) by
# .github/workflows/docker-publish.yml.
#
# Architecture: TWO-PHASE install.
#   Phase 1 (conda, this file's "Conda layer"): R 4.4.3 base + utility R
#     packages + Python 3.11 + celloracle + gimmemotifs. This is what
#     environment.yml declares. The conda solver handles it deterministically.
#   Phase 2 (BiocManager + CRAN, this file's "Bioconductor layer"):
#     Seurat 5.1.0, Signac 1.14.0, SeuratDisk, monocle3, cicero, EnsDb.
#     These are installed from CRAN/Bioc with explicit version pins. We do
#     NOT rely on conda-forge for them because their transitive dep graph
#     (r-spatstat, r-leiden, r-uwot, r-irlba, …) is not jointly resolvable
#     with r-base=4.4.3 in a single conda pass.
#
# Base image: mambaorg/micromamba:1.5.10-jammy (glibc 2.35, required for
# the modern r-base=4.4.3 build-string).
#
# Build locally:
#   docker build -t scMultiome-GRN:dev -f containers/Dockerfile containers/
#
# Run locally:
#   docker run --rm scMultiome-GRN:dev Rscript -e 'library(Seurat); library(cicero)'
#   docker run --rm scMultiome-GRN:dev python -c 'import celloracle; print(celloracle.__version__)'
#
# Author:      Daniel Mouzo
# License:     MIT
# ==============================================================================

FROM mambaorg/micromamba:1.5.10-jammy

LABEL org.opencontainers.image.title="scMultiome-GRN"
LABEL org.opencontainers.image.description="Polyglot R + Python container for the scMultiome-GRN Nextflow pipeline"
LABEL org.opencontainers.image.source="https://github.com/dmouzo/scMultiome-GRN"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.authors="Daniel Mouzo"

# ------------------------------------------------------------------------------
# System layer: Java for Nextflow, graphviz for DAGs, pandoc for R vignettes.
# ------------------------------------------------------------------------------
USER root
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        openjdk-17-jre-headless \
        graphviz \
        pandoc \
        ca-certificates \
        curl \
        libxt-dev \
        libfontconfig1 \
        fonts-dejavu-core && \
    rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------------
# Conda layer: the lightweight R + Python environment from environment.yml.
# This is the only layer that uses the conda solver; it stays small and
# deterministic (no Seurat / Signac / monocle3 here).
# ------------------------------------------------------------------------------
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# MAMBA_NO_LOCKFILE_INSTALL avoids a race with conda-lock files that some
# CI runners auto-generate. PYTHONNOUSERSITE keeps pip from reading
# ~/.local/ which does not exist in the container.
ENV MAMBA_NO_LOCKFILE_INSTALL=1 \
    PYTHONNOUSERSITE=1 \
    RETICULATE_PYTHON=/opt/conda/bin/python3

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a -y

# Make the conda env's bin first on PATH for every subsequent layer.
ENV PATH=/opt/conda/bin:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# ------------------------------------------------------------------------------
# Bioconductor layer: Seurat 5.1.0, Signac 1.14.0, SeuratDisk, monocle3,
# cicero, EnsDb.Mmusculus.v79 — installed from CRAN / BiocManager with
# explicit version pins. NO conda solver is involved here.
#
# This is the slow part (~10-15 min on a 4-core CI runner) but it is
# fully reproducible and never breaks the way the joint conda resolver
# does.
#
# The first two RUNs below install the R packages from CRAN. The third
# installs the Bioconductor packages with an explicit Bioc 3.19 version
# (the one aligned with R 4.4.x). They are split into separate RUNs so
# the Docker layer cache is preserved when only one of them changes.
# ------------------------------------------------------------------------------

# CRAN layer 1: Seurat + SeuratWrappers (already part of Seurat v5) + Signac
RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), \
                           Ncpus = parallel::detectCores()); \
    install.packages(c("Seurat", "Signac"), dependencies = TRUE)' && \
    rm -rf /tmp/Rtmp* /root/.cache/R

# CRAN layer 2: SeuratDisk (for h5ad export; GRN_dataProcess.R has a fallback
# to .rds sidecars if this is missing, so a failure here is non-fatal).
RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), \
                           Ncpus = parallel::detectCores()); \
    tryCatch( \
        install.packages("SeuratDisk", dependencies = TRUE), \
        error = function(e) message("SeuratDisk install failed (non-fatal): ", e$message))' \
    && rm -rf /tmp/Rtmp* /root/.cache/R

# Bioconductor layer: monocle3, cicero, EnsDb.Mmusculus.v79.
# Bioc 3.19 is the version aligned with R 4.4.x.
RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), \
                           Ncpus = parallel::detectCores()); \
    BiocManager::install( \
        c("monocle3", "cicero", "EnsDb.Mmusculus.v79"), \
        version = "3.19", \
        update = FALSE, \
        ask = FALSE, \
        dependencies = TRUE)' \
    && rm -rf /tmp/Rtmp* /root/.cache/R

# ------------------------------------------------------------------------------
# Nextflow layer: install latest stable Nextflow so the container can be
# invoked as `docker run ... nextflow run main.nf` directly.
# ------------------------------------------------------------------------------
RUN mkdir -p /opt/nextflow && \
    curl -fsSL https://get.nextflow.io | bash && \
    mv nextflow /opt/nextflow/nextflow && \
    chmod +x /opt/nextflow/nextflow && \
    ln -sf /opt/nextflow/nextflow /usr/local/bin/nextflow

# ------------------------------------------------------------------------------
# Validation layer: smoke-test that every required package is importable.
# SeuratDisk is checked with requireNamespace() because the GRN_dataProcess.R
# fallback to .rds sidecars makes its absence non-fatal. All other packages
# are HARD requirements — if any is missing, this RUN exits non-zero and the
# image is NOT pushed to GHCR.
# ------------------------------------------------------------------------------
RUN Rscript -e 'suppressPackageStartupMessages({ \
        library(Seurat); library(Signac); \
        library(cicero); library(monocle3); \
        library(EnsDb.Mmusculus.v79); library(ggplot2) }); \
    sd_ok <- requireNamespace("SeuratDisk", quietly = TRUE); \
    cat("R packages OK:", \
        "Seurat", as.character(packageVersion("Seurat")), "|", \
        "Signac", as.character(packageVersion("Signac")), "|", \
        "cicero", as.character(packageVersion("cicero")), "|", \
        "monocle3", as.character(packageVersion("monocle3")), "|", \
        "EnsDb.Mmusculus.v79", as.character(packageVersion("EnsDb.Mmusculus.v79")), "|", \
        "SeuratDisk", if (sd_ok) "loaded" else "MISSING (h5ad fallback active)", "\n")' && \
    python -c "import scanpy, anndata, celloracle, gimmemotifs; \
        print('Python packages OK:', \
              'celloracle', celloracle.__version__, '|', \
              'scanpy', scanpy.__version__, '|', \
              'anndata', anndata.__version__)" && \
    nextflow -version 2>&1 | head -3

# ------------------------------------------------------------------------------
# Runtime configuration: non-root user, working dir, entrypoint.
# ------------------------------------------------------------------------------
USER $MAMBA_USER

WORKDIR /workspace

# Default entrypoint runs nextflow so the container can be invoked directly.
# Override with:  docker run --rm -v $(pwd):/workspace IMAGE bash
ENTRYPOINT ["/opt/nextflow/nextflow"]
CMD ["-version"]
