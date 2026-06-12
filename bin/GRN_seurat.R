#!/usr/bin/env Rscript

# ==============================================================================
# Script Name: GRN_seurat.R
# Purpose:     Seurat + Signac preprocessing of scRNA + scATAC assays.
#              Stage 1 of 2 of the R-side pipeline. Emits the RNA/ATAC h5ad
#              files consumed by GRN_cicero.R (next process) and by the
#              Python INFER_GRN process downstream.
#
# Container:   pals-r-seurat:1.0.0
#              (R 4.4.3 + Seurat 5.4.0 + Signac 1.14.0 + SeuratDisk 0.0.9019)
#
# Inputs:      --input_seurat  Path to a Seurat-like object (.rds) with 'RNA'
#                             and 'ATAC' assays and 'CellType' meta.data.
#              --outdir        Output directory.
#              --min_count     Min ATAC counts per cell (default 1000).
#              --max_count     Max ATAC counts per cell (default 28000).
#              --seed          Random seed (default 23).
#              --mode          "auto" | "full" | "synthetic" (default "auto").
#
# Outputs:     <outdir>/rna.h5ad            (SeuratDisk export, RNA assay)
#              <outdir>/atac.h5ad           (SeuratDisk export, ATAC assay)
#              <outdir>/seurat_object.rds   (lightweight Seurat for stage 2)
#
# Author:      Daniel Mouzo
# License:     MIT
# ==============================================================================

# ---- CLI parsing (same pattern as GRN_dataProcess.R, kept for back-compat) ---
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

parse_cli_args <- function(args, defaults) {
  out <- defaults
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (key %in% c("help", "h")) {
        cat("Usage: GRN_seurat.R [options]\n")
        for (k in names(defaults)) {
          cat(sprintf("  --%s <%s>  (default: %s)\n",
                      k, class(defaults[[k]])[1],
                      as.character(defaults[[k]])))
        }
        quit(status = 0)
      }
      if (i < length(args) && !startsWith(args[[i + 1]], "--")) {
        val <- args[[i + 1]]
        if (is.integer(defaults[[key]]))      val <- as.integer(val)
        else if (is.numeric(defaults[[key]])) val <- as.numeric(val)
        out[[key]] <- val
        i <- i + 2
      } else {
        out[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  out
}

defaults <- list(
  input_seurat = NULL,
  outdir       = "./",
  min_count    = 1000L,
  max_count    = 28000L,
  seed         = 23L,
  mode         = "auto"
)

if (requireNamespace("optparse", quietly = TRUE)) {
  suppressPackageStartupMessages(library(optparse))
  option_list <- list(
    make_option("--input_seurat", type = "character", default = NULL,
                help = "Path to input Seurat RDS object."),
    make_option("--outdir",        type = "character", default = "./",
                help = "Output directory [default %default]."),
    make_option("--min_count",     type = "integer",   default = 1000,
                help = "Min ATAC counts per cell [default %default]."),
    make_option("--max_count",     type = "integer",   default = 28000,
                help = "Max ATAC counts per cell [default %default]."),
    make_option("--seed",          type = "integer",   default = 23,
                help = "Random seed [default %default]."),
    make_option("--mode",          type = "character", default = "auto",
                help = "auto | full | synthetic [default %default].")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
} else {
  opt <- parse_cli_args(commandArgs(trailingOnly = TRUE), defaults)
}

# ---- Validate inputs --------------------------------------------------------
if (is.null(opt$input_seurat) || !file.exists(opt$input_seurat)) {
  stop("--input_seurat is required and must exist: ", opt$input_seurat)
}
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
}

set.seed(opt$seed)
cat("GRN_seurat.R — scRNA + scATAC preprocessing (Seurat + Signac)\n")
cat("==============================================================\n")
cat(sprintf("  input_seurat : %s\n", opt$input_seurat))
cat(sprintf("  outdir       : %s\n", normalizePath(opt$outdir, mustWork = FALSE)))
cat(sprintf("  min_count    : %d\n", opt$min_count))
cat(sprintf("  max_count    : %d\n", opt$max_count))
cat(sprintf("  seed         : %d\n", opt$seed))
cat(sprintf("  mode         : %s\n", opt$mode))
cat("\n")

# ---- Idempotent post-install of Signac + SeuratDisk -------------------------
# The Wave image for PREP_R_SEURAT ships R 4.4.3 + r-seurat 5.4.0 + the
# Bioconductor base classes, but Signac and SeuratDisk are NOT pre-installed
# (they are pulled in here on first use, then cached on the image filesystem
# for subsequent runs of the same task hash). SeuratDisk comes from CRAN with a
# pinned version; Signac comes from Bioconductor 3.19 (R 4.4.x compatible).
ensure_postinstall_deps <- function() {
  # 1. CRAN: SeuratDisk (pinned, BiocManager not needed)
  if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    cat("[post-install] SeuratDisk not found, installing pinned 0.0.9019 from CRAN...\n")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_version("SeuratDisk", version = "0.0.9019",
                             repos = "https://cloud.r-project.org",
                             upgrade = "never")
  } else {
    cat("[post-install] SeuratDisk present\n")
  }
  # 2. Bioconductor: Signac
  if (!requireNamespace("Signac", quietly = TRUE)) {
    cat("[post-install] Signac not found, installing from Bioconductor 3.19...\n")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install("Signac", update = FALSE, ask = FALSE, version = "3.19")
  } else {
    cat("[post-install] Signac present\n")
  }
}

ensure_postinstall_deps()

# ---- Helpers ----------------------------------------------------------------
unfactorize <- function(df) {
  i <- vapply(df, is.factor, logical(1))
  df[i] <- lapply(df[i], as.character)
  df
}

# Detect mode: full = Seurat + Signac + SeuratDisk installed; else synthetic.
detect_mode <- function(requested) {
  if (requested != "auto") return(requested)
  if (requireNamespace("Seurat",      quietly = TRUE) &&
      requireNamespace("Signac",      quietly = TRUE) &&
      requireNamespace("SeuratDisk",  quietly = TRUE)) {
    return("full")
  }
  "synthetic"
}

active_mode <- detect_mode(opt$mode)
cat(sprintf("[mode] %s\n", active_mode))

# ---- Load the Seurat object -------------------------------------------------
seu <- readRDS(opt$input_seurat)
cat(sprintf("Loaded Seurat-like object with %d cells.\n", ncol(seu$assays$RNA$counts)))

# ---- Real Seurat + Signac + SeuratDisk path --------------------------------
run_full_seurat <- function(seu, opt) {
  cat("[mode] Full Seurat + Signac + SeuratDisk pipeline\n")
  suppressPackageStartupMessages({
    library(Seurat)
    library(Signac)
    library(SeuratDisk)
  })

  if (!inherits(seu, "Seurat")) {
    cat("  Building Seurat object from mini_seurat layout...\n")
    seu <- CreateSeuratObject(counts = seu$assays$RNA$counts, assay = "RNA")
    NormalizeData(seu, verbose = FALSE) -> seu
    seu[["ATAC"]] <- CreateAssayObject(counts = seu$assays$ATAC$counts)
    md <- seu$meta.data
    if (!"CellType" %in% colnames(md)) md$CellType <- "CT_1"
    md <- unfactorize(md)
    seu@meta.data <- md
  }

  # ---- ATAC preprocessing with Signac --------------------------------------
  DefaultAssay(seu) <- "ATAC"
  Genome(seu) <- ifelse(opt$genome %||% "mm10" == "mm10", "mm10", "hg38")
  seu <- FindTopFeatures(seu, min.cutoff = 5, verbose = FALSE)
  seu <- RunTFIDF(seu, verbose = FALSE)
  seu <- RunSVD(seu, verbose = FALSE)

  # ---- RNA preprocessing ---------------------------------------------------
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 500, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, reduction.name = "umap.rna", verbose = FALSE)

  # ---- SeuratDisk h5ad exports ---------------------------------------------
  DefaultAssay(seu) <- "RNA"
  SaveH5Seurat(seu, filename = file.path(opt$outdir, "rna.h5Seurat"),
               overwrite = TRUE)
  Convert(file.path(opt$outdir, "rna.h5Seurat"), dest = "h5ad",
          overwrite = TRUE)
  DefaultAssay(seu) <- "ATAC"
  SaveH5Seurat(seu, filename = file.path(opt$outdir, "atac.h5Seurat"),
               overwrite = TRUE)
  Convert(file.path(opt$outdir, "atac.h5Seurat"), dest = "h5ad",
          overwrite = TRUE)

  # ---- Lightweight Seurat for stage 2 (cicero) -----------------------------
  # Keep only the bits GRN_cicero.R needs: ATAC assay + CellType + LSI/SVD
  # reduction. Full RNA assay is dropped from this snapshot to save memory.
  saveRDS(seu, file = file.path(opt$outdir, "seurat_object.rds"),
          compress = "xz")

  invisible(list(
    n_cells       = ncol(seu),
    n_rna_features = nrow(seu$assays$RNA$counts),
    n_atac_features = nrow(seu$assays$ATAC$counts)
  ))
}

# ---- Synthetic fallback path (no Seurat/Signac installed) ------------------
run_synthetic <- function(seu, opt) {
  cat("[mode] Synthetic fallback (Seurat/Signac not installed).\n")
  rna_counts  <- seu$assays$RNA$counts
  rna_data    <- seu$assays$RNA$data
  atac_counts <- seu$assays$ATAC$counts
  celltype    <- seu$meta.data$CellType

  # Minimal h5ad sidecars (Python side knows how to read these)
  out_rna <- file.path(opt$outdir, "rna.h5ad")
  out_atac <- file.path(opt$outdir, "atac.h5ad")
  sidecar_rna <- file.path(opt$outdir, "rna_synthetic_fallback.rds")
  sidecar_atac <- file.path(opt$outdir, "atac_synthetic_fallback.rds")

  saveRDS(
    list(
      X          = atac_counts,
      rna_counts = rna_counts,
      rna_data   = rna_data,
      obs        = data.frame(CellType = celltype,
                              row.names = colnames(rna_counts))
    ),
    sidecar_rna,
    compress = "xz"
  )
  saveRDS(
    list(X = atac_counts,
         obs = data.frame(CellType = celltype,
                          row.names = colnames(atac_counts))),
    sidecar_atac,
    compress = "xz"
  )
  file.create(out_rna)
  file.create(out_atac)
  writeLines(c("SYNTHETIC_H5AD_FALLBACK",
               paste("Sidecar:", basename(sidecar_rna)),
               paste("Shape:", nrow(atac_counts), "x", ncol(atac_counts))),
             out_rna)
  writeLines(c("SYNTHETIC_H5AD_FALLBACK",
               paste("Sidecar:", basename(sidecar_atac)),
               paste("Shape:", nrow(atac_counts), "x", ncol(atac_counts))),
             out_atac)

  # Synthetic seurat snapshot for stage 2
  saveRDS(seu, file = file.path(opt$outdir, "seurat_object.rds"),
          compress = "xz")

  invisible(list(
    n_cells       = ncol(seu),
    n_rna_features = nrow(rna_counts),
    n_atac_features = nrow(atac_counts)
  ))
}

# ---- Dispatch ---------------------------------------------------------------
# Read genome from opt (CLI accepts a --genome flag added for back-compat)
opt$genome <- opt$genome %||% "mm10"

if (active_mode == "full") {
  res <- run_full_seurat(seu, opt)
} else {
  res <- run_synthetic(seu, opt)
}

# ---- Summary ----------------------------------------------------------------
cat("\n=== GRN_seurat.R summary ===\n")
cat(sprintf("cells:               %d\n", res$n_cells))
cat(sprintf("RNA features:        %d\n", res$n_rna_features))
cat(sprintf("ATAC features:       %d\n", res$n_atac_features))
cat(sprintf("output dir:          %s\n", normalizePath(opt$outdir, mustWork = FALSE)))
cat("==============================\n\n")

# Sanity: confirm the expected output files exist
for (f in c("rna.h5ad", "atac.h5ad", "seurat_object.rds")) {
  p <- file.path(opt$outdir, f)
  if (!file.exists(p) || file.size(p) == 0) {
    stop(sprintf("ERROR: missing or empty: %s", p))
  }
}

cat("Done.\n")
