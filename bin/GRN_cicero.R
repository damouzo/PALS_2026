#!/usr/bin/env Rscript

# ==============================================================================
# Script Name: GRN_cicero.R
# Purpose:     Cicero co-accessibility computation over the ATAC assay.
#              Stage 2 of 2 of the R-side pipeline. Consumes the lightweight
#              Seurat object emitted by GRN_seurat.R and emits the deterministic
#              all_peaks.csv + cicero_connections.csv tables consumed by the
#              Python INFER_GRN process.
#
# Container:   pals-r-bioc:1.0.0
#              (R 4.4.3 + Bioconductor 3.19 + cicero + monocle3 +
#               EnsDb.Mmusculus.v79)
#
# Inputs:      --seurat_rds    Path to the lightweight Seurat object emitted
#                             by GRN_seurat.R (file named seurat_object.rds).
#              --outdir        Output directory (same as the one passed to
#                             GRN_seurat.R; cicero tables are written here).
#              --min_count     Min ATAC counts per cell (default 1000).
#              --max_count     Max ATAC counts per cell (default 28000).
#              --seed          Random seed (default 23).
#              --mode          "auto" | "full" | "synthetic" (default "auto").
#
# Outputs:     <outdir>/all_peaks.csv           (one column, chr-prefixed peaks)
#              <outdir>/cicero_connections.csv  (Peak1, Peak2, coaccess)
#
# Author:      Daniel Mouzo
# License:     MIT
# ==============================================================================

# ---- CLI parsing ------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

parse_cli_args <- function(args, defaults) {
  out <- defaults
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (key %in% c("help", "h")) {
        cat("Usage: GRN_cicero.R [options]\n")
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
  seurat_rds = NULL,
  outdir     = "./",
  min_count  = 1000L,
  max_count  = 28000L,
  seed       = 23L,
  mode       = "auto"
)

if (requireNamespace("optparse", quietly = TRUE)) {
  suppressPackageStartupMessages(library(optparse))
  option_list <- list(
    make_option("--seurat_rds", type = "character", default = NULL,
                help = "Path to the lightweight Seurat RDS from GRN_seurat.R."),
    make_option("--outdir",     type = "character", default = "./",
                help = "Output directory [default %default]."),
    make_option("--min_count",  type = "integer",   default = 1000,
                help = "Min ATAC counts per cell [default %default]."),
    make_option("--max_count",  type = "integer",   default = 28000,
                help = "Max ATAC counts per cell [default %default]."),
    make_option("--seed",       type = "integer",   default = 23,
                help = "Random seed [default %default]."),
    make_option("--mode",       type = "character", default = "auto",
                help = "auto | full | synthetic [default %default].")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
} else {
  opt <- parse_cli_args(commandArgs(trailingOnly = TRUE), defaults)
}

# ---- Validate inputs --------------------------------------------------------
if (is.null(opt$seurat_rds) || !file.exists(opt$seurat_rds)) {
  stop("--seurat_rds is required and must exist: ", opt$seurat_rds)
}
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
}

set.seed(opt$seed)
cat("GRN_cicero.R — Cicero co-accessibility\n")
cat("=======================================\n")
cat(sprintf("  seurat_rds : %s\n", opt$seurat_rds))
cat(sprintf("  outdir     : %s\n", normalizePath(opt$outdir, mustWork = FALSE)))
cat(sprintf("  min_count  : %d\n", opt$min_count))
cat(sprintf("  max_count  : %d\n", opt$max_count))
cat(sprintf("  seed       : %d\n", opt$seed))
cat(sprintf("  mode       : %s\n", opt$mode))
cat("\n")

# ---- Helpers ----------------------------------------------------------------
detect_mode <- function(requested) {
  if (requested != "auto") return(requested)
  if (requireNamespace("cicero",   quietly = TRUE) &&
      requireNamespace("monocle3", quietly = TRUE)) {
    return("full")
  }
  "synthetic"
}

active_mode <- detect_mode(opt$mode)
cat(sprintf("[mode] %s\n", active_mode))

# ---- Real cicero + monocle3 path -------------------------------------------
run_full_cicero <- function(seu, opt) {
  cat("[mode] Full cicero + monocle3 pipeline\n")
  suppressPackageStartupMessages({
    library(cicero)
    library(monocle3)
    # Use the M. musculus v79 Ensembl annotation shipped in the bioc image.
    # The mm10 assembly is the default for the microglia demo dataset.
    library(EnsDb.Mmusculus.v79)
  })

  DefaultAssay(seu) <- "ATAC"

  atac_cds <- cicero::as.cell_data_set(seu)
  if ("UMAP_ATAC" %in% names(reducedDims(atac_cds))) {
    reducedDims(atac_cds)$UMAP <- reducedDims(atac_cds)$UMAP_ATAC
    reducedDims(atac_cds)$UMAP_ATAC <- NULL
  }

  atac_cds <- cicero::detect_genes(atac_cds)
  atac_cds <- atac_cds[Matrix::rowSums(SummarizedExperiment::exprs(atac_cds)) != 0, ]
  atac_cds <- atac_cds[, Matrix::colSums(SummarizedExperiment::exprs(atac_cds)) >= opt$min_count]
  atac_cds <- atac_cds[, Matrix::colSums(SummarizedExperiment::exprs(atac_cds)) <= opt$max_count]
  atac_cds <- cicero::estimate_size_factors(atac_cds)
  atac_cds <- cicero::preprocess_cds(atac_cds, method = "LSI")
  atac_cds <- cicero::reduce_dimension(atac_cds, reduction_method = "UMAP",
                                       preprocess_method = "LSI")
  umap_coords <- reducedDims(atac_cds)$UMAP
  cicero_cds <- cicero::make_cicero_cds(atac_cds,
                                        reduced_coordinates = umap_coords)

  ref <- seqlengths(EnsDb.Mmusculus.v79)
  chromosome_length <- data.frame(V1 = names(ref), V2 = ref)
  rownames(chromosome_length) <- seq_len(nrow(chromosome_length))
  chromosome_length <- chromosome_length[nchar(chromosome_length$V1) <= 2, ]

  conns <- cicero::run_cicero(cicero_cds, chromosome_length)

  all_peaks <- rownames(SummarizedExperiment::exprs(atac_cds))
  write.csv(all_peaks,
            file = file.path(opt$outdir, "all_peaks.csv"),
            row.names = FALSE, quote = FALSE)
  write.csv(conns,
            file = file.path(opt$outdir, "cicero_connections.csv"),
            row.names = FALSE, quote = FALSE)

  # ---- Memory hygiene ------------------------------------------------------
  rm(seu, atac_cds, cicero_cds, conns); gc(verbose = FALSE)

  invisible(list(n_peaks = length(all_peaks),
                 n_connections = nrow(conns)))
}

# ---- Synthetic fallback path (no cicero installed) -------------------------
run_synthetic_cicero <- function(seu, opt) {
  cat("[mode] Synthetic fallback pipeline (cicero not installed).\n")
  atac_counts <- seu$assays$ATAC$counts

  n_peaks <- nrow(atac_counts)
  set.seed(opt$seed)
  anchors <- sample.int(n_peaks, size = min(50, n_peaks))
  partners <- unlist(lapply(anchors, function(a) {
    cand <- setdiff(seq_len(n_peaks), a)
    sample(cand, size = min(10, length(cand)))
  }))
  rows <- rep(anchors, each = 10)[seq_along(partners)]
  cols <- partners
  keep <- rows < cols
  rows <- rows[keep]; cols <- cols[keep]
  pair_key <- paste(pmin(rows, cols), pmax(rows, cols), sep = "-")
  keep <- !duplicated(pair_key)
  rows <- rows[keep]; cols <- cols[keep]

  if (length(rows) > 0) {
    ac1 <- atac_counts[rows, , drop = FALSE]
    ac2 <- atac_counts[cols, , drop = FALSE]
    cor_mat <- sapply(seq_along(rows), function(i) {
      v1 <- as.numeric(ac1[i, ]); v2 <- as.numeric(ac2[i, ])
      if (sd(v1) == 0 || sd(v2) == 0) return(0)
      suppressWarnings(cor(v1, v2, method = "spearman"))
    })
    coaccess <- round(pmax(0, cor_mat), 4)
  } else {
    coaccess <- numeric(0)
  }

  conns <- data.frame(
    Peak1    = rownames(atac_counts)[rows],
    Peak2    = rownames(atac_counts)[cols],
    coaccess = coaccess,
    stringsAsFactors = FALSE
  )
  conns <- conns[conns$coaccess > 0, ]

  all_peaks <- rownames(atac_counts)
  write.csv(all_peaks,
            file = file.path(opt$outdir, "all_peaks.csv"),
            row.names = FALSE, quote = FALSE)
  write.csv(conns,
            file = file.path(opt$outdir, "cicero_connections.csv"),
            row.names = FALSE, quote = FALSE)

  invisible(list(n_peaks = length(all_peaks),
                 n_connections = nrow(conns)))
}

# ---- Dispatch ---------------------------------------------------------------
seu <- readRDS(opt$seurat_rds)

if (active_mode == "full") {
  res <- run_full_cicero(seu, opt)
} else {
  res <- run_synthetic_cicero(seu, opt)
}

# ---- Summary ----------------------------------------------------------------
cat("\n=== GRN_cicero.R summary ===\n")
cat(sprintf("peaks exported:        %d\n", res$n_peaks))
cat(sprintf("connections exported:  %d\n", res$n_connections))
cat(sprintf("output dir:            %s\n", normalizePath(opt$outdir, mustWork = FALSE)))
cat("==================================\n\n")

# Sanity
for (f in c("all_peaks.csv", "cicero_connections.csv")) {
  p <- file.path(opt$outdir, f)
  if (!file.exists(p) || file.size(p) == 0) {
    stop(sprintf("ERROR: missing or empty: %s", p))
  }
}

cat("Done.\n")
