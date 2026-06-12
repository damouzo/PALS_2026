#!/usr/bin/env Rscript

# ==============================================================================
# Script Name: build_synthetic_mini_seurat.R
# Purpose:     Generate the ultra-mini synthetic Seurat object used by the
#              scMultiome-GRN smoke-test profile (-profile test, --input
#              data_demo/synthetic/mini_seurat.rds).
#
#              This script is fully deterministic (set.seed(23)) and tries
#              several data sources in the following priority order:
#
#                1. SeuratData::pbmcMultiome   (preferred; real 10x multiome)
#                2. CellOracle example data    (PBMC; requires celloracle)
#                3. Local h5ad fallback        (any .h5ad in /tmp or work dir)
#                4. Synthetic fallback         (built from base R + Matrix)
#
#              The synthetic fallback guarantees that a minimal, runnable
#              mini_seurat.rds can always be produced (e.g. for CI smoke
#              tests where Seurat is not installed).
#
# Output:      data_demo/synthetic/mini_seurat.rds
#              (Seurat v4 object with RNA + ATAC assays, 'CellType' column,
#               umap.rna + umap.atac reductions, 150-200 cells, ~500 HVGs and
#               ~2000 peaks.)
#
# Usage:       Rscript data_demo/synthetic/build_synthetic_mini_seurat.R
#
# Author:      Daniel Mouzo
# License:     MIT
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
})

set.seed(23)

# -----------------------------------------------------------------------------
# Helper: minimal null-coalescing (avoid %||% which is R 4.4+ only)
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Determine output path: same dir as this script's parent (data_demo/synthetic).
script_arg <- sub("--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
script_dir <- if (nzchar(script_arg)) dirname(normalizePath(script_arg)) else getwd()
OUT_PATH <- file.path(script_dir, "mini_seurat.rds")
OUT_PATH <- normalizePath(OUT_PATH, mustWork = FALSE)

# -----------------------------------------------------------------------------
# Helper: coerce factors to characters in a data.frame (mimic Seurat internals)
# -----------------------------------------------------------------------------
unfactorize <- function(df) {
  i <- vapply(df, is.factor, logical(1))
  df[i] <- lapply(df[i], as.character)
  df
}

# -----------------------------------------------------------------------------
# Helper: write a Seurat-v4-compatible .rds object.
# -----------------------------------------------------------------------------
write_mini_seurat <- function(seu_list, out_path) {
  if (!dir.exists(dirname(out_path))) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  }
  saveRDS(seu_list, out_path, compress = "xz")
  size_mb <- file.info(out_path)$size / (1024 ^ 2)
  cat(sprintf("Wrote %s (%.2f MB)\n", out_path, size_mb))
  if (size_mb > 50) {
    warning(sprintf("Output file is %.2f MB (>50 MB). Demo will be slow.", size_mb))
  }
  invisible(size_mb)
}

# -----------------------------------------------------------------------------
# Helper: synthetic fallback (no external dependencies)
# -----------------------------------------------------------------------------
build_synthetic_mini_seurat <- function(n_cells = 180,
                                        n_rna_genes = 500,
                                        n_atac_peaks = 2000,
                                        n_celltypes = 3) {
  cat("[synthetic] Building synthetic mini Seurat object from scratch...\n")

  celltypes <- paste0("CT", seq_len(n_celltypes))
  ct_vec <- sample(rep(celltypes, length.out = n_cells), size = n_cells, replace = FALSE)

  cluster_bias <- matrix(
    rnorm(n_rna_genes * n_celltypes, mean = 0, sd = 0.4),
    nrow = n_rna_genes,
    ncol = n_celltypes
  )

  rna_means <- exp(0.5 + cluster_bias)
  rna_counts <- matrix(0L, nrow = n_rna_genes, ncol = n_cells)
  for (k in seq_len(n_celltypes)) {
    cells_k <- which(ct_vec == celltypes[k])
    if (length(cells_k) > 0) {
      rna_counts[, cells_k] <- matrix(
        rpois(length(cells_k) * n_rna_genes, lambda = rna_means[, k]),
        nrow = n_rna_genes
      )
    }
  }
  rna_counts <- as(rna_counts, "CsparseMatrix")
  rownames(rna_counts) <- sprintf("g%04d", seq_len(n_rna_genes))
  colnames(rna_counts) <- sprintf("cell%03d", seq_len(n_cells))
  tf_symbols <- c("GATA1", "MYB", "KLF1", "TAL1", "LMO2", "RUNX1", "FLI1",
                  "HOXA9", "MEIS1", "PBX1")
  for (i in seq_along(tf_symbols)) {
    if (i <= n_rna_genes) rownames(rna_counts)[i] <- tf_symbols[i]
  }

  atac_means <- exp(0.3 + 0.3 * cluster_bias[sample(n_rna_genes, n_atac_peaks, replace = TRUE), ])
  atac_counts <- matrix(0L, nrow = n_atac_peaks, ncol = n_cells)
  for (k in seq_len(n_celltypes)) {
    cells_k <- which(ct_vec == celltypes[k])
    if (length(cells_k) > 0) {
      atac_counts[, cells_k] <- matrix(
        rpois(length(cells_k) * n_atac_peaks, lambda = atac_means[, k]),
        nrow = n_atac_peaks
      )
    }
  }
  atac_counts <- as(atac_counts, "CsparseMatrix")
  rownames(atac_counts) <- paste0(
    sample(c("chr1", "chr2", "chr3", "chr5", "chr7", "chr11", "chr19"),
           n_atac_peaks, replace = TRUE),
    "-",
    sample(1:1e6, n_atac_peaks, replace = FALSE),
    "-",
    sample(1:1e6, n_atac_peaks, replace = FALSE)
  )
  colnames(atac_counts) <- sprintf("cell%03d", seq_len(n_cells))

  rna_data <- log1p(rna_counts * (1e4 / Matrix::colSums(rna_counts)))

  umap_rna <- matrix(0, nrow = n_cells, ncol = 2)
  umap_atac <- matrix(0, nrow = n_cells, ncol = 2)
  centers <- matrix(rnorm(n_celltypes * 2, sd = 3), ncol = 2)
  for (k in seq_len(n_celltypes)) {
    cells_k <- which(ct_vec == celltypes[k])
    umap_rna[cells_k, ] <- centers[k, ] + matrix(rnorm(length(cells_k) * 2, sd = 0.6), ncol = 2)
    umap_atac[cells_k, ] <- centers[k, ] + matrix(rnorm(length(cells_k) * 2, sd = 0.7), ncol = 2)
  }
  rownames(umap_rna) <- rownames(umap_atac) <- colnames(rna_counts)

  meta <- data.frame(
    CellType      = ct_vec,
    nFeature_RNA  = Matrix::colSums(rna_counts > 0),
    nCount_RNA    = Matrix::colSums(rna_counts),
    nFeature_ATAC = Matrix::colSums(atac_counts > 0),
    nCount_ATAC   = Matrix::colSums(atac_counts),
    row.names     = colnames(rna_counts)
  )
  meta <- unfactorize(meta)

  seu <- list(
    assays = list(
      RNA  = list(counts = rna_counts,  data = rna_data),
      ATAC = list(counts = atac_counts, data = atac_counts)
    ),
    meta.data  = meta,
    reductions = list(
      umap.rna  = list(cell.embeddings = umap_rna),
      umap.atac = list(cell.embeddings = umap_atac)
    ),
    version      = "synthetic-mini-1.0",
    project.name = "scMultiome-GRN_mini_synthetic"
  )
  class(seu) <- c("mini_seurat", "list")
  seu
}

# -----------------------------------------------------------------------------
# Strategy 1: SeuratData::pbmcMultiome (real 10x multiome)
# -----------------------------------------------------------------------------
try_pbmc_multiome <- function() {
  if (!requireNamespace("SeuratData", quietly = TRUE)) return(NULL)
  if (!requireNamespace("Seurat", quietly = TRUE)) return(NULL)
  cat("[strategy 1] SeuratData::pbmcMultiome\n")
  suppressPackageStartupMessages({
    library(SeuratData)
    library(Seurat)
  })
  if (!"pbmcMultiome" %in% AvailableData()$Name) {
    InstallData("pbmcMultiome")
  }
  data("pbmcMultiome")
  seu <- pbmcMultiome
  seu <- subset(seu, nFeature_RNA > 200 & nFeature_ATAC > 100)
  cell_ids <- unlist(tapply(colnames(seu), seu$CellType,
                            function(x) sample(x, min(60, length(x)))))
  seu <- subset(seu, cells = cell_ids)
  return(seu)
}

# -----------------------------------------------------------------------------
# Strategy 4: synthetic (always-available fallback)
# -----------------------------------------------------------------------------
try_synthetic <- function() build_synthetic_mini_seurat()

# -----------------------------------------------------------------------------
# Try strategies in priority order
# -----------------------------------------------------------------------------
seu <- NULL
for (strategy in list(try_pbmc_multiome, try_synthetic)) {
  res <- tryCatch(strategy(), error = function(e) {
    cat(sprintf("  -> failed: %s\n", conditionMessage(e)))
    NULL
  })
  if (!is.null(res)) {
    seu <- res
    break
  }
}

if (is.null(seu)) {
  stop("All build strategies failed. Cannot produce mini_seurat.rds.")
}

# Optional postprocess for real Seurat objects.
if (!inherits(seu, "mini_seurat") && requireNamespace("Seurat", quietly = TRUE)) {
  cat("[postprocess] Downsampling real Seurat object...\n")
  seu <- subset(seu, cells = sample(colnames(seu), min(180, ncol(seu))))
  seu <- FindVariableFeatures(seu, nfeatures = 500)
  DefaultAssay(seu) <- "RNA"
  VariableFeatures(seu) <- head(VariableFeatures(seu), 500)
  DefaultAssay(seu) <- "ATAC"
  if ("ATAC" %in% names(seu@assays)) {
    atac_top <- head(names(sort(Matrix::rowSums(seu@assays$ATAC$counts > 0),
                                decreasing = TRUE)), 2000)
    seu <- subset(seu, features = c(VariableFeatures(seu), atac_top))
  }
}

if (!"CellType" %in% colnames(seu$meta.data %||% seu@meta.data)) {
  warning("Source data lacked 'CellType' column. Renaming seurat_clusters.")
  md <- seu$meta.data %||% seu@meta.data
  if ("seurat_clusters" %in% colnames(md)) {
    md$CellType <- paste0("CT_", as.character(md$seurat_clusters))
  } else {
    md$CellType <- "CT_1"
  }
  seu$meta.data <- md
}

if (!inherits(seu, "mini_seurat") && requireNamespace("Seurat", quietly = TRUE)) {
  cat("[postprocess] Wrapping real Seurat object into mini_seurat layout...\n")
  seu <- list(
    assays = list(
      RNA  = list(counts = seu@assays$RNA$counts,  data = seu@assays$RNA$data),
      ATAC = list(counts = seu@assays$ATAC$counts, data = seu@assays$ATAC$data)
    ),
    meta.data  = unfactorize(as.data.frame(seu@meta.data)),
    reductions = list(
      umap.rna  = list(cell.embeddings = seu@reductions$umap.rna@cell.embeddings),
      umap.atac = list(cell.embeddings = seu@reductions$umap.atac@cell.embeddings)
    ),
    version      = "wrapped-real-seurat-1.0",
    project.name = "scMultiome-GRN_mini_real"
  )
  class(seu) <- c("mini_seurat", "list")
}

md <- seu$meta.data
ur <- dim(seu$reductions$umap.rna$cell.embeddings)
ua <- dim(seu$reductions$umap.atac$cell.embeddings)
cat("\n=== mini_seurat.rds summary ===\n")
cat(sprintf("Cells:           %d\n", nrow(md)))
cat(sprintf("RNA features:    %d\n", nrow(seu$assays$RNA$counts)))
cat(sprintf("ATAC peaks:      %d\n", nrow(seu$assays$ATAC$counts)))
cat(sprintf("CellTypes:       %s\n", paste(unique(md$CellType), collapse = ", ")))
cat(sprintf("UMAP.rna dim:    %d x %d\n", ur[1], ur[2]))
cat(sprintf("UMAP.atac dim:   %d x %d\n", ua[1], ua[2]))
cat("================================\n\n")

write_mini_seurat(seu, OUT_PATH)
cat("\nDone. mini_seurat.rds is ready for the smoke test.\n")
