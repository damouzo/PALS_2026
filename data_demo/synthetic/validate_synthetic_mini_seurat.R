#!/usr/bin/env Rscript

# ==============================================================================
# Script Name: validate_synthetic_mini_seurat.R
# Purpose:     Sanity-check data_demo/synthetic/mini_seurat.rds.
#              Verifies all the fields the downstream pipeline expects.
# ==============================================================================

suppressPackageStartupMessages({ library(Matrix) })

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

script_arg <- sub("--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
script_dir <- if (nzchar(script_arg)) dirname(normalizePath(script_arg)) else getwd()
RDS_PATH <- file.path(script_dir, "mini_seurat.rds")
RDS_PATH <- normalizePath(RDS_PATH, mustWork = FALSE)

if (!file.exists(RDS_PATH)) {
  stop("mini_seurat.rds not found at ", RDS_PATH)
}

cat("Loading", RDS_PATH, "...\n")
seu <- readRDS(RDS_PATH)
size_mb <- file.info(RDS_PATH)$size / (1024 ^ 2)
cat(sprintf("  size: %.2f MB\n", size_mb))

required <- list(
  "assays$RNA$counts"  = function(s) nrow(s$assays$RNA$counts)  > 0 && ncol(s$assays$RNA$counts)  > 0,
  "assays$RNA$data"    = function(s) !is.null(s$assays$RNA$data),
  "assays$ATAC$counts" = function(s) nrow(s$assays$ATAC$counts) > 0 && ncol(s$assays$ATAC$counts) > 0,
  "meta.data$CellType" = function(s) "CellType" %in% colnames(s$meta.data),
  "reductions$umap.rna"  = function(s) !is.null(s$reductions$umap.rna$cell.embeddings),
  "reductions$umap.atac" = function(s) !is.null(s$reductions$umap.atac$cell.embeddings)
)

stop_flag <- FALSE
for (fld in names(required)) {
  ok <- tryCatch(required[[fld]](seu), error = function(e) FALSE)
  cat(sprintf("  %-30s %s\n", fld, if (ok) "OK" else "FAIL"))
  if (!ok) stop_flag <- TRUE
}

n_cells   <- ncol(seu$assays$RNA$counts)
n_rna     <- nrow(seu$assays$RNA$counts)
n_atac    <- nrow(seu$assays$ATAC$counts)
n_cts     <- length(unique(seu$meta.data$CellType))

cat("\n--- Numeric checks ---\n")
cat(sprintf("  cells in [150, 200]:   %s (%d)\n", if (n_cells >= 150 && n_cells <= 200) "OK" else "FAIL", n_cells))
cat(sprintf("  RNA features == 500:   %s (%d)\n", if (n_rna == 500) "OK" else "FAIL", n_rna))
cat(sprintf("  ATAC peaks == 2000:    %s (%d)\n", if (n_atac == 2000) "OK" else "FAIL", n_atac))
cat(sprintf("  celltypes >= 2:        %s (%d)\n", if (n_cts >= 2) "OK" else "FAIL", n_cts))
cat(sprintf("  file size < 50 MB:     %s (%.2f MB)\n", if (size_mb < 50) "OK" else "FAIL", size_mb))

cat("\n--- Cell type distribution ---\n")
print(table(seu$meta.data$CellType))

tf_present <- intersect(c("GATA1", "MYB", "KLF1"), rownames(seu$assays$RNA$counts))
cat(sprintf("\n--- TF targets detected: %s ---\n", paste(tf_present, collapse = ", ")))

if (stop_flag) stop("One or more required fields missing.")
cat("\nAll checks passed.\n")
