#!/usr/bin/env python3

# ==============================================================================
# Script Name: GRN_analysis.py
# Purpose:     Build a base Gene Regulatory Network with CellOracle, score
#              network edges, and simulate in silico transcription-factor
#              knockouts. All plots are saved to --outdir.
#
#              The script supports two execution modes:
#                - "full" : real CellOracle + GimmeMotifs pipeline.
#                - "synthetic" : deterministic fallback that reproduces the
#                  exact same output schema (base_GRN.parquet, oracle.h5,
#                  links.h5, plots/*.png) for CI smoke tests and environments
#                  where CellOracle is not installed.
#
# Inputs:      --peaks         all_peaks.csv from GRN_dataProcess.R
#              --connections   cicero_connections.csv from GRN_dataProcess.R
#              --rna_h5ad      rna.h5ad (real) or rna_synthetic_fallback.rds
#              --atac_h5ad     atac.h5ad (real) or atac_synthetic_fallback.rds
#              --percentile    Network link pruning percentile (WOW factor)
#              --motif_db      GimmeMotifs motif database
#              --genome        Reference genome assembly (hg38, mm10, ...)
#              --clustering    Cell-type column in the AnnData object
#              --umap_key      UMAP embedding key in the AnnData object
#              --clust2study   Cluster to use for in silico perturbation
#              --ko_targets    Comma-separated TF symbols to simulate KO
#              --outdir        Output directory
#              --mode          "auto", "full" or "synthetic"
#
# Outputs:     <outdir>/base_GRN.parquet
#              <outdir>/oracle.celloracle.oracle
#              <outdir>/links.celloracle.links
#              <outdir>/network_scores.xlsx
#              <outdir>/plots/vector_field_KO_<TF>.png     [per KO target]
#              <outdir>/plots/network_ranked_score.png
#              <outdir>/plots/score_comparison.png
#
# Author:      Daniel Mouzo
# License:     MIT
# ==============================================================================

from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

# Matplotlib must be non-interactive for CI / no-display environments.
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Silence noisy FutureWarnings from scanpy / anndata.
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("grn_analysis")


# ==============================================================================
# Resilient serialization helpers
# ==============================================================================
_PARQUET_ENGINE = None
try:
    import pyarrow  # noqa: F401
    _PARQUET_ENGINE = "pyarrow"
except ImportError:
    try:
        import fastparquet  # noqa: F401
        _PARQUET_ENGINE = "fastparquet"
    except ImportError:
        pass


def save_dataframe(df: pd.DataFrame, path: Path, fmt: Optional[str] = None) -> Path:
    """Save a DataFrame to parquet (preferred) or CSV (fallback)."""
    fmt = fmt or ("parquet" if _PARQUET_ENGINE else "csv")
    if fmt == "parquet" and _PARQUET_ENGINE:
        try:
            df.to_parquet(path, engine=_PARQUET_ENGINE)
            return path
        except Exception as e:
            log.warning("parquet write failed (%s); falling back to CSV", e)
            fmt = "csv"
    csv_path = path.with_suffix(".csv") if path.suffix == ".parquet" else path
    df.to_csv(csv_path, index=True)
    return csv_path


# ==============================================================================
# CLI parsing
# ==============================================================================
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="GRN_analysis.py",
        description=(
            "Infer a base Gene Regulatory Network with CellOracle, score the "
            "network and simulate in silico TF knockouts."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--peaks",        required=True, type=str,
                   help="all_peaks.csv from GRN_dataProcess.R")
    p.add_argument("--connections",  required=True, type=str,
                   help="cicero_connections.csv from GRN_dataProcess.R")
    p.add_argument("--rna_h5ad",     required=True, type=str,
                   help="rna.h5ad (real) or rna_synthetic_fallback.rds")
    p.add_argument("--atac_h5ad",    required=True, type=str,
                   help="atac.h5ad (real) or atac_synthetic_fallback.rds")
    p.add_argument("--percentile",   required=True, type=int, default=90,
                   help="Network link pruning percentile (the -resume demo param).")
    p.add_argument("--motif_db",     type=str, default="gimme.vertebrate.v5.0",
                   help="GimmeMotifs motif database identifier.")
    p.add_argument("--genome",       type=str, default="hg38",
                   help="Reference genome assembly.")
    p.add_argument("--clustering",   type=str, default="CellType",
                   help="Cell-type column in the AnnData object.")
    p.add_argument("--umap_key",     type=str, default="umap.rna",
                   help="UMAP embedding key in the AnnData object.")
    p.add_argument("--clust2study",  type=str, default="CT1",
                   help="Cluster to use for in silico perturbation.")
    p.add_argument("--ko_targets",   type=str, default="GATA1,MYB,KLF1",
                   help="Comma-separated TF symbols to simulate as knockouts.")
    p.add_argument("--outdir",       type=str, default="./",
                   help="Output directory.")
    p.add_argument("--mode",         type=str, default="auto",
                   choices=["auto", "full", "synthetic"],
                   help="Execution mode (auto detects from environment).")
    p.add_argument("--p_threshold",  type=float, default=0.001,
                   help="P-value threshold for link filtering.")
    p.add_argument("--n_top_edges",  type=int, default=2000,
                   help="Maximum number of edges to keep after filtering.")
    p.add_argument("--seed",         type=int, default=23,
                   help="Random seed for reproducibility.")
    return p


# ==============================================================================
# Data loading helpers (support both real h5ad and synthetic fallback sidecars)
# ==============================================================================
def _read_synthetic_sidecar(path: str) -> dict:
    """Read a synthetic_fallback.rds file written by GRN_dataProcess.R."""
    import pyreadr  # lightweight; pulls rds via pandas
    try:
        result = pyreadr.read_r(path)
        # pyreadr returns a dict of DataFrames; the synthetic file has a
        # single key. Just return the raw pyreadr dict.
        return result
    except ImportError:
        # Fallback: try to parse with anndata (pyreadr not installed).
        log.warning("pyreadr not installed; using pandas fallback for .rds sidecar.")
        return {}


def load_atac_anndata(atac_path: str, atac_synthetic_sidecar: Optional[str]) -> "anndata.AnnData":
    """Load the ATAC AnnData. Supports real .h5ad and synthetic fallback."""
    if atac_path.endswith(".h5ad") and os.path.getsize(atac_path) > 1024:
        import anndata as ad
        adata = sc.read_h5ad(atac_path) if False else ad.read_h5ad(atac_path)
        return adata
    # Synthetic fallback path: read the sidecar.
    if atac_synthetic_sidecar and os.path.exists(atac_synthetic_sidecar):
        try:
            import anndata as ad
            import scipy.sparse as sp
            result = _read_synthetic_sidecar(atac_synthetic_sidecar)
            if not result:
                return _synthetic_anndata()
            df = list(result.values())[0]
            X = sp.csr_matrix(df.values.astype(float))
            obs = pd.DataFrame(index=df.index.astype(str))
            var = pd.DataFrame(index=df.columns.astype(str))
            return ad.AnnData(X=X, obs=obs, var=var)
        except Exception as e:
            log.warning("Failed to read sidecar %s: %s", atac_synthetic_sidecar, e)
    return _synthetic_anndata()


def load_rna_anndata(rna_path: str, rna_synthetic_sidecar: Optional[str]) -> "anndata.AnnData":
    """Load the RNA AnnData. Supports real .h5ad and synthetic fallback."""
    if rna_path.endswith(".h5ad") and os.path.getsize(rna_path) > 1024:
        import anndata as ad
        return ad.read_h5ad(rna_path)
    if rna_synthetic_sidecar and os.path.exists(rna_synthetic_sidecar):
        try:
            import anndata as ad
            import scipy.sparse as sp
            result = _read_synthetic_sidecar(rna_synthetic_sidecar)
            if not result:
                return _synthetic_anndata()
            # The R sidecar has 'rna_counts' and 'rna_data' as separate keys.
            counts_df = result.get("rna_counts")
            if counts_df is not None:
                X = sp.csr_matrix(counts_df.values.astype(float))
                var = pd.DataFrame(index=counts_df.columns.astype(str))
                obs = pd.DataFrame(index=counts_df.index.astype(str))
                adata = ad.AnnData(X=X, obs=obs, var=var)
                if "rna_data" in result:
                    norm_df = result["rna_data"]
                    adata.X = sp.csr_matrix(norm_df.values.astype(float))
                return adata
        except Exception as e:
            log.warning("Failed to read sidecar %s: %s", rna_synthetic_sidecar, e)
    return _synthetic_anndata()


def _synthetic_anndata() -> "anndata.AnnData":
    """Last-resort: a tiny AnnData with random expression for the UMAP plot."""
    import anndata as ad
    import scipy.sparse as sp
    n_cells, n_genes = 50, 30
    X = sp.csr_matrix(np.random.poisson(1.0, size=(n_cells, n_genes)))
    obs = pd.DataFrame({"CellType": ["CT1"] * 25 + ["CT2"] * 25},
                       index=[f"cell{i:03d}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"g{i:04d}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


# ==============================================================================
# Synthetic fallback pipeline (no CellOracle installed)
# ==============================================================================
def run_synthetic_pipeline(args: argparse.Namespace,
                           rna_counts: pd.DataFrame,
                           atac_counts: pd.DataFrame,
                           celltype: pd.Series) -> None:
    """Reproduce the schema of the full CellOracle pipeline using only
    numpy/pandas/matplotlib. Used for CI smoke tests and offline demos.
    """
    log.info("Mode: synthetic fallback (CellOracle not installed).")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plots = outdir / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    np.random.seed(args.seed)

    # ---- 1. Base GRN (TF -> target weight matrix) ----
    tf_symbols = ["GATA1", "MYB", "KLF1", "TAL1", "LMO2", "RUNX1", "FLI1", "HOXA9", "MEIS1", "PBX1"]
    tfs = [t for t in tf_symbols if t in rna_counts.index]
    targets = rna_counts.index.difference(tfs).tolist()
    log.info("TFs detected: %d (%s); targets: %d", len(tfs), ", ".join(tfs), len(targets))

    n_tfs, n_targets = len(tfs), min(300, len(targets))
    selected_targets = np.random.choice(targets, size=n_targets, replace=False)
    base_grn = pd.DataFrame(
        np.random.lognormal(mean=0.0, sigma=1.0, size=(n_tfs, n_targets)).astype("float32"),
        index=tfs,
        columns=selected_targets,
    )
    base_grn.index.name = "source"
    save_dataframe(base_grn, outdir / "base_GRN.parquet")
    log.info("Wrote base_GRN.parquet (%d TFs x %d targets)", n_tfs, n_targets)

    # ---- 2. Per-cluster network scores ----
    clusters = sorted(celltype.unique())
    n_edges = min(args.n_top_edges, n_tfs * n_targets)
    rows = []
    for cl in clusters:
        # Per-cluster score: simulated cluster-specific edge weights
        cluster_weights = base_grn.values * np.random.lognormal(0, 0.3, base_grn.shape)
        score_df = pd.DataFrame(cluster_weights, index=base_grn.index,
                                columns=base_grn.columns)
        # Compute centrality-like metric: row sum (outgoing strength)
        centrality = score_df.sum(axis=1).sort_values(ascending=False)
        for rank, (tf, val) in enumerate(centrality.items(), start=1):
            rows.append({"cluster": cl, "rank": rank, "TF": tf,
                         "degree_centrality_out": val})
    score_df_full = pd.DataFrame(rows)

    # Excel: one sheet per cluster
    try:
        with pd.ExcelWriter(outdir / "network_scores.xlsx",
                            engine="openpyxl") as writer:
            for cl in clusters:
                cl_df = score_df_full[score_df_full["cluster"] == cl]
                cl_df.to_excel(writer, sheet_name=str(cl)[:31], index=False)
        log.info("Wrote network_scores.xlsx")
    except ImportError:
        # openpyxl missing; fall back to plain CSV
        score_df_full.to_csv(outdir / "network_scores.csv", index=False)
        log.info("Wrote network_scores.csv (openpyxl missing)")

    # ---- 3. Link filtering (the WOW-factor path) ----
    threshold = np.percentile(score_df_full["degree_centrality_out"].values,
                              args.percentile)
    filtered = score_df_full[score_df_full["degree_centrality_out"] >= threshold].copy()
    log.info("Percentile %d -> threshold %.4f -> %d links kept",
             args.percentile, threshold, len(filtered))

    # ---- 4. UMAP-based static visualizations ----
    # Compute a 2D embedding from RNA counts (PCA -> 2D) for plotting.
    from sklearn.decomposition import PCA
    dense = rna_counts.values if hasattr(rna_counts, "values") else np.asarray(rna_counts)
    dense_log = np.log1p(dense)
    n_comp = min(2, min(dense_log.shape) - 1)
    if n_comp >= 2:
        umap = PCA(n_components=2).fit_transform(dense_log.T)
    else:
        umap = np.random.RandomState(args.seed).randn(dense_log.shape[1], 2)
    umap_df = pd.DataFrame(umap, columns=["UMAP1", "UMAP2"],
                           index=rna_counts.columns)
    umap_df["CellType"] = celltype.reindex(umap_df.index).values

    # Plot 1: ranked network score
    fig, ax = plt.subplots(figsize=(7, 4.5))
    top = filtered.sort_values("degree_centrality_out", ascending=False).head(30)
    if len(top) > 0:
        ax.barh(top["TF"].astype(str), top["degree_centrality_out"], color="steelblue")
        ax.invert_yaxis()
    ax.set_xlabel("degree_centrality_out")
    ax.set_title(f"Network ranked score (cluster={args.clust2study}, "
                 f"percentile={args.percentile})")
    fig.tight_layout()
    fig.savefig(plots / "network_ranked_score.png", dpi=150)
    plt.close(fig)
    log.info("Wrote plots/network_ranked_score.png")

    # Plot 2: score comparison between two clusters
    if len(clusters) >= 2:
        c1, c2 = clusters[0], clusters[1]
        s1 = score_df_full[score_df_full["cluster"] == c1].set_index("TF")["degree_centrality_out"]
        s2 = score_df_full[score_df_full["cluster"] == c2].set_index("TF")["degree_centrality_out"]
        common = s1.index.intersection(s2.index)
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.scatter(s1[common], s2[common], s=18, alpha=0.6, c="darkgreen")
        ax.set_xlabel(c1)
        ax.set_ylabel(c2)
        ax.set_title(f"Score comparison: {c1} vs {c2} (pctl={args.percentile})")
        fig.tight_layout()
        fig.savefig(plots / "score_comparison.png", dpi=150)
        plt.close(fig)
        log.info("Wrote plots/score_comparison.png")

    # ---- 5. ⭐ In silico perturbation: vector-field plots per KO target ----
    ko_list = [t.strip() for t in args.ko_targets.split(",") if t.strip()]
    for ko in ko_list:
        fig, ax = plt.subplots(figsize=(7, 6))
        # Plot baseline cells
        for cl, color in zip(clusters, ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]):
            mask = umap_df["CellType"] == cl
            ax.scatter(umap_df.loc[mask, "UMAP1"],
                       umap_df.loc[mask, "UMAP2"],
                       s=12, alpha=0.5, label=cl, c=color)
        # Vector field: synthetic gradient pointing "away" from the KO
        # cluster centroid. The sign of the shift is derived from the
        # base GRN weight of `ko` across all targets.
        if ko in base_grn.index:
            shift = base_grn.loc[ko].values
            shift = (shift - shift.mean()) / (shift.std() + 1e-9)
            # Pick two principal directions
            direction = np.array([shift[0], shift[len(shift) // 2]])
            direction /= np.linalg.norm(direction) + 1e-9
        else:
            direction = np.array([1.0, 0.0])
        xx, yy = np.meshgrid(np.linspace(umap_df["UMAP1"].min(), umap_df["UMAP1"].max(), 12),
                             np.linspace(umap_df["UMAP2"].min(), umap_df["UMAP2"].max(), 12))
        u = direction[0] * np.ones_like(xx)
        v = direction[1] * np.ones_like(yy)
        ax.quiver(xx, yy, u, v, alpha=0.5, color="black",
                  scale=15, width=0.004, headwidth=3)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title(f"In silico KO of {ko}: predicted cell-fate shift vector field")
        ax.legend(loc="best", fontsize=8, frameon=False)
        fig.tight_layout()
        out_path = plots / f"vector_field_KO_{ko}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        log.info("Wrote plots/%s", out_path.name)

    # ---- 6. Serialize "oracle" and "links" objects as parquet placeholders ----
    # The downstream tooling expects these names. We emit parquet files
    # that summarize the run.
    oracle_summary = pd.DataFrame({
        "cluster": clusters,
        "n_cells": [int((celltype == cl).sum()) for cl in clusters],
    })
    save_dataframe(oracle_summary, outdir / "oracle.celloracle.oracle")
    save_dataframe(filtered, outdir / "links.celloracle.links")
    log.info("Wrote oracle.celloracle.oracle and links.celloracle.links")


# ==============================================================================
# Real CellOracle pipeline (used when celloracle is importable)
# ==============================================================================
def run_full_pipeline(args: argparse.Namespace,
                      rna_adata, atac_adata) -> None:
    """End-to-end CellOracle + GimmeMotifs pipeline. Mirrors the canonical
    CellOracle tutorial (Kamimoto et al., 2023).
    """
    log.info("Mode: full CellOracle + GimmeMotifs pipeline.")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plots = outdir / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    import scanpy as sc
    import celloracle as co
    from celloracle import motif_analysis as ma
    from gimmemotifs.motif import read_motifs, MotifConfig

    # ---- 1. Load data and prepare layers ----
    atac_adata.layers["counts"] = atac_adata.raw.X.copy() if atac_adata.raw is not None else atac_adata.X.copy()
    atac_adata.layers["norm"]   = atac_adata.X.copy()
    rna_adata.layers["counts"]  = rna_adata.raw.X.copy() if rna_adata.raw is not None else rna_adata.X.copy()
    rna_adata.layers["norm"]    = rna_adata.X.copy()

    peaks = pd.read_csv(args.peaks, header=0).iloc[:, 0].tolist()
    peaks = np.array(
        ["chr" + str(p).replace("-", "_") for p in peaks],
        dtype=object,
    )
    conns = pd.read_csv(args.connections, index_col=0)

    def _norm(v: object) -> str:
        v = str(v).replace("-", "_")
        return v if v.startswith("chr") else "chr" + v
    conns[["Peak1", "Peak2"]] = conns[["Peak1", "Peak2"]].applymap(_norm)

    # ---- 2. TSS annotation ----
    tss = ma.get_tss_info(peak_str_list=peaks, ref_genome=args.genome)
    integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss,
                                                   cicero_connections=conns)
    peaks_filtered = integrated[integrated.coaccess >= 1][
        ["peak_id", "gene_short_name"]
    ].reset_index(drop=True)
    peaks_filtered.to_csv(outdir / "processed_peak_file.csv", index=False)

    # ---- 3. Motif scanning ----
    peaks_checked = ma.check_peak_format(peaks_filtered, args.genome)
    tfi = ma.TFinfo(peak_data_frame=peaks_checked, ref_genome=args.genome)
    config = MotifConfig()
    motif_dir = config.get_motif_dir()
    motif_path = os.path.join(motif_dir, args.motif_db + ".pfm")
    motifs = read_motifs(motif_path)
    tfi.scan(motifs=motifs, verbose=True)
    tfi.to_hdf5(file_path=str(outdir / "tfinfo.celloracle.tfinfo"))

    # ---- 4. Base GRN ----
    tfi.reset_filtering()
    tfi.filter_motifs_by_score(threshold=10)
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
    base_grn = tfi.to_dataframe()
    save_dataframe(base_grn, outdir / "base_GRN.parquet")

    # ---- 5. Oracle object ----
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=rna_adata,
        cluster_column_name=args.clustering,
        embedding_name=args.umap_key,
    )
    oracle.import_TF_data(TF_info_matrix=base_grn)

    oracle.perform_PCA()
    n_comps = min(np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_) > 0.002)))[0][0], 50)
    k = int(0.025 * oracle.adata.shape[0])
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True,
                          b_sight=k * 8, b_maxl=k * 4, n_jobs=4)
    oracle.to_hdf5(str(outdir / "oracle.celloracle.oracle"))

    # ---- 6. GRN calculation ----
    links = oracle.get_links(cluster_name_for_GRN_unit=args.clustering,
                              alpha=10, verbose_level=10)
    links.to_hdf5(file_path=str(outdir / "links.celloracle.links"))

    # ---- 7. Network filtering (the WOW-factor path) ----
    links.filter_links(p=args.p_threshold, weight="coef_abs",
                       threshold_number=args.n_top_edges)
    links.filter_links(percentile=args.percentile)
    links.get_network_score()
    try:
        links.merged_score.to_excel(str(outdir / "network_scores.xlsx"))
    except ImportError:
        save_dataframe(links.merged_score.reset_index(),
                       outdir / "network_scores.xlsx")

    # ---- 8. Static network plots ----
    links.plot_scores_as_rank(cluster=args.clust2study, n_gene=30,
                              save=str(plots / "network_ranked_score"))
    clusters_available = list(links.cluster)
    if len(clusters_available) >= 2:
        links.plot_score_comparison_2D(
            value="degree_centrality_all",
            cluster1=clusters_available[0], cluster2=clusters_available[1],
            percentile=args.percentile,
            save=str(plots / "score_comparison"),
        )
    plt.close("all")

    # ---- 9. ⭐ In silico perturbation: vector fields ----
    ko_list = [t.strip() for t in args.ko_targets.split(",") if t.strip()]
    for ko in ko_list:
        try:
            oracle.simulate_shift(ko_target=ko, n_propagation=3)
            links.plot_simulation_flow(ko_target=ko,
                                       save=str(plots / f"vector_field_KO_{ko}"))
        except Exception as e:
            log.warning("KO simulation for %s failed: %s", ko, e)
        plt.close("all")

    log.info("Full CellOracle pipeline completed.")


# ==============================================================================
# Mode dispatch
# ==============================================================================
def detect_mode(args: argparse.Namespace) -> str:
    if args.mode != "auto":
        return args.mode
    try:
        import celloracle  # noqa: F401
        return "full"
    except ImportError:
        return "synthetic"


def main() -> None:
    args = build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    np.random.seed(args.seed)

    log.info("GRN_analysis.py — CellOracle GRN inference + perturbation")
    log.info("=========================================================")
    for k, v in vars(args).items():
        log.info("  %-15s : %s", k, v)
    log.info("")

    # Sidecar detection (R fallback path)
    rna_sidecar = (args.rna_h5ad if args.rna_h5ad.endswith(".rds")
                   else str(Path(args.rna_h5ad).with_name(
                       Path(args.rna_h5ad).stem + "_synthetic_fallback.rds")))
    atac_sidecar = (args.atac_h5ad if args.atac_h5ad.endswith(".rds")
                    else str(Path(args.atac_h5ad).with_name(
                        Path(args.atac_h5ad).stem + "_synthetic_fallback.rds")))

    mode = detect_mode(args)
    log.info("Selected mode: %s", mode)

    if mode == "full":
        atac_adata = load_atac_anndata(args.atac_h5ad, atac_sidecar)
        rna_adata  = load_rna_anndata(args.rna_h5ad, rna_sidecar)
        run_full_pipeline(args, rna_adata, atac_adata)
    else:
        # Synthetic fallback loads the R-sidecar via the original mini_seurat
        # object. To avoid round-tripping through R, we read the sidecar
        # parquet/rds directly.
        try:
            rna = pd.read_csv  # placeholder to test imports
        except Exception:
            pass
        # Reconstruct minimal dataframes from the synthetic sidecars.
        rna_counts, atac_counts, celltype = _load_synthetic_inputs(rna_sidecar, atac_sidecar)
        run_synthetic_pipeline(args, rna_counts, atac_counts, celltype)

    log.info("Done. Outputs in: %s", outdir.resolve())


def _load_synthetic_inputs(rna_sidecar: str, atac_sidecar: str
                            ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """Load the R synthetic sidecars and project them into the DataFrames
    the synthetic Python pipeline expects.
    """
    try:
        import pyreadr
        rna_obj = pyreadr.read_r(rna_sidecar)
        atac_obj = pyreadr.read_r(atac_sidecar)
        # The R sidecar contains: rna_counts (cells x genes), rna_data,
        # obs (with CellType column). For atac, only X is meaningful.
        rna_counts_df = list(rna_obj["rna_counts"].values)
        # We need a (genes x cells) DataFrame, transpose
        rna_counts = pd.DataFrame(rna_obj["rna_counts"]).T
        rna_counts.index.name = "gene"
        atac_df = pd.DataFrame(atac_obj.get("X",
                                            list(atac_obj.values())[0]))
        atac_counts = atac_df  # (peaks x cells)
        # Cell type
        obs = rna_obj.get("obs")
        if obs is None and "obs" not in rna_obj:
            # Last resort
            celltype = pd.Series(["CT1"] * rna_counts.shape[1],
                                  index=rna_counts.columns)
        else:
            obs_df = list(rna_obj.values())[0] if obs is None else obs
            celltype = pd.Series(
                obs_df["CellType"].values if "CellType" in obs_df.columns
                else ["CT1"] * len(obs_df),
                index=rna_counts.columns,
            )
        log.info("Loaded synthetic sidecars: RNA=%s ATAC=%s",
                 rna_counts.shape, atac_counts.shape)
        return rna_counts, atac_counts, celltype
    except ImportError:
        log.warning("pyreadr not installed; building fully synthetic inputs.")
        return _pure_synthetic_inputs()


def _pure_synthetic_inputs() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    np.random.seed(23)
    n_cells, n_genes, n_peaks = 180, 500, 2000
    tf_symbols = ["GATA1", "MYB", "KLF1", "TAL1", "LMO2", "RUNX1", "FLI1", "HOXA9", "MEIS1", "PBX1"]
    gene_names = tf_symbols + [f"g{i:04d}" for i in range(n_genes - len(tf_symbols))]
    peak_names = [f"chr{np.random.choice([1,2,3,5,7,11,19])}-{a}-{b}"
                  for a, b in zip(np.random.randint(1, 1e6, n_peaks),
                                  np.random.randint(1, 1e6, n_peaks))]
    cell_names = [f"cell{i:03d}" for i in range(n_cells)]
    rna_counts = pd.DataFrame(
        np.random.poisson(1.0, size=(n_genes, n_cells)),
        index=gene_names, columns=cell_names,
    )
    atac_counts = pd.DataFrame(
        np.random.poisson(0.5, size=(n_peaks, n_cells)),
        index=peak_names, columns=cell_names,
    )
    celltype = pd.Series(
        np.random.choice(["CT1", "CT2", "CT3"], size=n_cells),
        index=cell_names,
    )
    return rna_counts, atac_counts, celltype


if __name__ == "__main__":
    main()
