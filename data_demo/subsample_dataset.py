#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# Script Name: subsample_dataset.py
# Purpose:     Create a lightweight, biologically meaningful test slice from
#              the 10x Genomics scMultiome Mouse Brain Alzheimer's AppNote
#              dataset, designed to feed the scMultiome-GRN Nextflow pipeline.
#
#              Uses the OFFICIAL Cell Ranger ARC 2.0.1 cluster labels
#              (graphclust) shipped in Multiome_..._analysis.tar.gz:
#                 - GEX cluster 21  = microglia (homeostatic + DAM-1)
#                 - GEX cluster 22  = DAM-like microglia (DAM-2)
#                 - ATAC cluster 7  = the cross-omics-equivalent microglia cluster
#              (Source: 10x Genomics AppNote, 10x_App-Note_Alzheimers_Letter.)
#
#              Default design: 3 CellTypes x 200 cells = 600 nuclei.
#                 - Microglia_Homeostatic  (WT 13.2m,  cluster 21)
#                 - Microglia_EarlyDAM     (TgCRND8 5.7m,  clusters 21 + 22)
#                 - Microglia_LateDAM      (TgCRND8 17.9m, cluster 22)
#              ATAC peaks are restricted to chr7 (Apoe), chr15 (Slc1a3) and
#              chr17 (Trem2) — the three regulatory loci highlighted in the
#              10x AppNote. Fragments are also filtered to those chromosomes
#              via tabix (with a pure-Python fallback).
#
#              The output is a Seurat v4-compatible .rds (via SeuratDisk)
#              AND an .h5ad (drop-in for the Python stack), plus a
#              filtered atac_fragments file and a run summary JSON.
#
# Inputs (downloaded by data_demo/download_data.sh):
#   data_demo/raw/Multiome_..._filtered_feature_bc_matrix.h5
#   data_demo/raw/Multiome_..._atac_fragments.tsv.gz
#   data_demo/raw/Multiome_..._atac_fragments.tsv.gz.tbi
#   data_demo/raw/Multiome_..._analysis.tar.gz  (10x cluster labels)
#   data_demo/raw/Multiome_..._web_summary.html
#   data_demo/raw/Multiome_..._summary.csv
#
# Outputs (in data_demo/processed/):
#   microglia_dam_demo.h5ad
#   microglia_dam_demo.rds
#   atac_fragments_chr7_chr15_chr17.tsv.gz
#   atac_fragments_chr7_chr15_chr17.tsv.gz.tbi  (if tabix available)
#   subsample_summary.json
#
# Usage:
#   python data_demo/subsample_dataset.py
#   python data_demo/subsample_dataset.py --help
#   python data_demo/subsample_dataset.py --cohorts my_cohorts.yaml
#   python data_demo/subsample_dataset.py --cells-per-cohort 100
#   python data_demo/subsample_dataset.py --no-rds --no-tabix
#
# Requirements (provided by containers/environment.yml):
#   scanpy, anndata, pandas, numpy, scipy, h5py
#   Optional: tabix (system) for fast fragment filtering
#   For .rds conversion: Rscript + SeuratDisk (run inside the container)
#
# Author:      Daniel Mouzo
# License:     MIT
# ==============================================================================

from __future__ import annotations

import argparse
import gzip
import io
import json
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd
    from scipy import sparse

# Heavy third-party imports are deferred to _import_heavy() so that
# `--help` works even on a host that does not have scanpy/anndata/numpy.

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("subsample_dataset")


# ==============================================================================
# Constants
# ==============================================================================
PREFIX = "Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote"

# Barcode suffix -> sample name (verified against summary.csv and analysis.tar.gz).
# WT timepoints: 2.5, 5.7, 13.2 months. TgCRND8: 2.5, 5.7, 17.9 months.
SUFFIX_TO_SAMPLE: Dict[int, str] = {
    3:  "AD_17p9_rep4",
    4:  "AD_17p9_rep5",
    5:  "AD_2p5_rep2",
    6:  "AD_2p5_rep3",
    7:  "AD_5p7_rep2",
    8:  "AD_5p7_rep6",
    9:  "WT_13p2_rep2",
    10: "WT_13p2_rep5",
    11: "WT_2p5_rep2",
    12: "WT_2p5_rep7",
    13: "WT_5p7_rep2",
    14: "WT_5p7_rep3",
}

# OFFICIAL 10x GEX graphclust cluster IDs (from analysis.tar.gz).
# Cluster 21 = microglia (11 of 12 canonical markers top-expressed here,
#   including Trem2, Tmem119, C1qa, Hexb, C1qb, C1qc, P2ry12, Cx3cr1, Csf1r,
#   Cd68, P2ry13). Cluster 22 = DAM-like, with P2ry12 specifically.
GEX_MICROGLIA_CLUSTER = 21
GEX_DAM_CLUSTER       = 22
ATAC_MICROGLIA_CLUSTER = 7  # cross-omics equivalent (820/846 GEX-mg barcodes)


# ==============================================================================
# Default cohort spec (used when --cohorts is not given)
# ==============================================================================
# Each cohort:
#   - label       : human-readable name printed in the log
#   - celltype    : final value of the CellType meta.data column
#   - suffixes    : list of barcode suffixes to draw from
#   - gex_clusters: official 10x GEX cluster IDs that are "on-target"
#   - n           : target number of nuclei
DEFAULT_COHORTS: List[Dict[str, Any]] = [
    {
        "label": "WT_13p2_microglia_homeostatic",
        "celltype": "Microglia_Homeostatic",
        "suffixes": [9, 10],
        "gex_clusters": [GEX_MICROGLIA_CLUSTER],
        "n": 200,
    },
    {
        "label": "AD_5p7_microglia_earlyDAM",
        "celltype": "Microglia_EarlyDAM",
        "suffixes": [7, 8],
        "gex_clusters": [GEX_MICROGLIA_CLUSTER, GEX_DAM_CLUSTER],
        "n": 200,
    },
    {
        "label": "AD_17p9_microglia_lateDAM",
        "celltype": "Microglia_LateDAM",
        "suffixes": [3, 4],
        "gex_clusters": [GEX_DAM_CLUSTER],
        "n": 200,
    },
]


# ==============================================================================
# Heavy-deps lazy loader
# ==============================================================================
_HEAVY: Dict[str, Any] = {}


def _import_heavy() -> Dict[str, Any]:
    if _HEAVY:
        return _HEAVY
    missing: List[str] = []
    for modname, pkg in [
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("sparse", "scipy"),
        ("sc", "scanpy"),
        ("anndata", "anndata"),
        ("h5py", "h5py"),
    ]:
        try:
            if modname in ("sparse",):
                from scipy import sparse as _sp
                _HEAVY[modname] = _sp
            elif modname in ("sc",):
                # scanpy is conventionally imported as `sc`; resolve the
                # real module name and store it under the alias key.
                _HEAVY[modname] = __import__("scanpy")
            elif modname in ("numpy", "pandas", "anndata", "h5py"):
                _HEAVY[modname] = __import__(modname)
            else:
                _HEAVY[modname] = __import__(modname)
        except ImportError:
            missing.append(pkg)
    if missing:
        log.error("Missing required Python packages: " + ", ".join(missing))
        log.error("Install them with:")
        log.error("    pip install " + " ".join(missing))
        log.error("Or run this script inside the Python GRN container:")
        log.error("    docker login ghcr.io                                    # one-time")
        log.error("    docker run --rm -u $(id -u):$(id -g) -v $(pwd):/workspace \\")
        log.error("        ghcr.io/damouzo/pals-python-grn:1.0.0 \\")
        log.error("        python /workspace/data_demo/subsample_dataset.py")
        sys.exit(1)
    return _HEAVY


# ==============================================================================
# CLI parsing
# ==============================================================================
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create a biologically meaningful microglia-only test "
                    "slice from the 10x Multiome Mouse Brain Alzheimer's dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--raw-dir", type=Path,
                   default=Path(__file__).parent / "raw",
                   help="Directory holding the raw 10x assets.")
    p.add_argument("--out-dir", type=Path,
                   default=Path(__file__).parent / "processed",
                   help="Directory for the processed outputs.")
    p.add_argument("--cohorts", type=Path, default=None,
                   help="Path to a YAML or JSON file overriding the default "
                        "cohort spec. Each entry must contain: label, celltype, "
                        "suffixes (list[int]), gex_clusters (list[int]), n (int). "
                        "If omitted, the built-in 3-cohort design is used.")
    p.add_argument("--cells-per-cohort", type=int, default=None,
                   help="Override the 'n' field of every cohort in the spec. "
                        "Useful for quick scaling (e.g. --cells-per-cohort 100).")
    p.add_argument("--target-chroms", nargs="+",
                   default=["chr7", "chr15", "chr17"],
                   help="Chromosomes to retain in the ATAC peak set and the "
                        "filtered fragments file. Defaults to chr7 (Apoe), "
                        "chr15 (Slc1a3, the distal peak highlighted in the 10x "
                        "AppNote) and chr17 (Trem2).")
    p.add_argument("--microglia-markers", nargs="+",
                   default=["C1qa", "Hexb", "Tmem119", "P2ry12"],
                   help="Canonical mouse microglia marker genes (used as a "
                        "QC report, not as a hard gate).")
    p.add_argument("--use-official-10x-clusters", action="store_true",
                   default=True,
                   help="[DEFAULT] Use the official Cell Ranger ARC graphclust "
                        "labels from analysis.tar.gz.")
    p.add_argument("--no-official-clusters", dest="use_official_10x_clusters",
                   action="store_false",
                   help="Disable the official 10x cluster labels and fall back "
                        "to a MicrogliaScore (mean log1p of the markers) gate. "
                        "Useful for testing the script on data without the "
                        "analysis.tar.gz available.")
    p.add_argument("--seed", type=int, default=23,
                   help="Random seed for deterministic subsampling.")
    p.add_argument("--no-rds", action="store_true",
                   help="Skip the Seurat .rds conversion (saves time if R is "
                        "not installed locally).")
    p.add_argument("--no-tabix", action="store_true",
                   help="Filter the ATAC fragments in pure Python instead of "
                        "delegating to the system `tabix` tool.")
    return p.parse_args()


# ==============================================================================
# Helpers
# ==============================================================================
def human_size(n: int) -> str:
    if n >= 1 << 30: return f"{n / (1 << 30):.2f} GB"
    if n >= 1 << 20: return f"{n / (1 << 20):.2f} MB"
    if n >= 1 << 10: return f"{n / (1 << 10):.2f} KB"
    return f"{n} B"


def banner(step: int, title: str) -> None:
    log.info("")
    log.info("=" * 72)
    log.info(f" STEP {step}: {title}")
    log.info("=" * 72)


def require_file(path: Path, what: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        log.error(f"Required {what} not found or empty: {path}")
        log.error("Run data_demo/download_data.sh first.")
        sys.exit(1)
    log.info(f"  {what:30s} {path}  ({human_size(path.stat().st_size)})")


def load_cohort_spec(args: argparse.Namespace) -> List[Dict[str, Any]]:
    """Resolve the active cohort spec from --cohorts / --cells-per-cohort."""
    spec: List[Dict[str, Any]]
    if args.cohorts is not None:
        log.info(f"Loading cohort spec from {args.cohorts}")
        text = args.cohorts.read_text()
        if args.cohorts.suffix.lower() in (".yaml", ".yml"):
            try:
                import yaml  # type: ignore
                spec = yaml.safe_load(text)
            except ImportError:
                log.error("PyYAML is required to load a YAML cohort spec. "
                          "Use a JSON file or `pip install pyyaml`.")
                sys.exit(1)
        else:
            spec = json.loads(text)
        log.info(f"  Loaded {len(spec)} cohort(s) from user spec.")
    else:
        spec = [dict(c) for c in DEFAULT_COHORTS]
        log.info(f"Using built-in default cohort spec ({len(spec)} cohorts).")
    # Apply --cells-per-cohort override uniformly.
    if args.cells_per_cohort is not None:
        for c in spec:
            c["n"] = args.cells_per_cohort
        log.info(f"  --cells-per-cohort override applied (n={args.cells_per_cohort}).")
    # Validate.
    required_keys = {"label", "celltype", "suffixes", "gex_clusters", "n"}
    for i, c in enumerate(spec):
        missing = required_keys - set(c.keys())
        if missing:
            log.error(f"Cohort {i} missing keys: {missing}")
            sys.exit(1)
    # Pretty-print.
    log.info("  Active cohort design:")
    log.info(f"    {'#':<3} {'label':<40} {'celltype':<28} {'n':>5}  suffixes / clusters")
    for i, c in enumerate(spec, 1):
        log.info(f"    {i:<3} {c['label']:<40} {c['celltype']:<28} "
                 f"{c['n']:>5}  {c['suffixes']} / {c['gex_clusters']}")
    return spec


# ==============================================================================
# Step 1: Load the aggregated 10x H5
# ==============================================================================
def load_h5(h5_path: Path):
    banner(1, "Loading aggregated 10x H5 matrix")
    sc = _import_heavy()["sc"]
    log.info(f"  Reading: {h5_path.name}")
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()
    if "feature_types" not in adata.var.columns:
        log.error("H5 lacks 'feature_types' column. Not a 10x multiome H5?")
        sys.exit(1)
    adata.var["modality"] = adata.var["modality"].fillna("Unknown") if "modality" in adata.var.columns else adata.var["feature_types"].map({
        "Gene Expression": "RNA",
        "Peaks":          "ATAC",
    }).fillna("Unknown")
    n_rna  = int((adata.var["modality"] == "RNA").sum())
    n_atac = int((adata.var["modality"] == "ATAC").sum())
    log.info(f"  Cells         : {adata.n_obs:,}")
    log.info(f"  Features      : {adata.n_vars:,}  "
             f"(RNA={n_rna:,}  ATAC={n_atac:,})")
    return adata


# ==============================================================================
# Step 2: Cohort / cluster assignment
# ==============================================================================
def annotate_barcodes(adata, spec: List[Dict[str, Any]],
                      use_official: bool) -> Dict[str, pd.DataFrame]:
    """Annotate each cell with its sample, cohort, and (optionally) 10x cluster.

    Returns a dict {barcode: row-dict} ready to attach as adata.obs.
    """
    np = _import_heavy()["np"]
    pd = _import_heavy()["pandas"]
    banner(2, "Assigning sample + cohort + (optional) 10x cluster labels")

    suffix = (
        adata.obs_names.to_series()
        .str.split("-").str[-1]
        .astype(int)
    )
    adata.obs["sample_suffix"] = suffix.values
    adata.obs["sample"] = adata.obs["sample_suffix"].map(SUFFIX_TO_SAMPLE)

    # Build the per-cohort candidate masks now (we still need clusters to
    # apply them; the official-cluster labels are added below).
    cohort_candidates: Dict[str, List[str]] = {}
    for c in spec:
        mask = adata.obs["sample_suffix"].isin(c["suffixes"])
        cohort_candidates[c["label"]] = adata.obs.index[mask].tolist()
        log.info(f"  {c['label']:40s}  sample suffixes={c['suffixes']}  "
                 f"raw cells={len(cohort_candidates[c['label']])}")
    return cohort_candidates


def load_10x_clusters(adata, analysis_tar: Path) -> Dict[str, int]:
    """Read the GEX graphclust clusters.csv from analysis.tar.gz without
    extracting the whole archive."""
    banner(2.5, "Loading OFFICIAL 10x GEX cluster labels (graphclust)")
    log.info(f"  Reading clusters.csv from {analysis_tar.name}")
    clusters: Dict[str, int] = {}
    n_target = adata.n_obs
    with tarfile.open(analysis_tar, "r:gz") as tf:
        member = None
        for m in tf.getmembers():
            if (m.name.endswith("clustering/gex/graphclust/clusters.csv")
                    and m.size < 5_000_000):
                member = m
                break
        if member is None:
            log.error("Could not find gex/graphclust/clusters.csv in tarball.")
            sys.exit(1)
        f = tf.extractfile(member)
        if f is None:
            log.error("Failed to extract clusters.csv.")
            sys.exit(1)
        # Skip header.
        next(f)
        for line in io.TextIOWrapper(f, encoding="utf-8"):
            line = line.rstrip("\n")
            if not line:
                continue
            barcode, cluster = line.split(",", 1)
            try:
                clusters[barcode] = int(cluster)
            except ValueError:
                continue
    log.info(f"  Loaded cluster assignments for {len(clusters):,} barcodes")
    adata.obs["gex_cluster_10x"] = (
        adata.obs_names.map(clusters).astype("Int64")
    )
    n_assigned = int(adata.obs["gex_cluster_10x"].notna().sum())
    log.info(f"  Barcodes in our H5 that got a cluster: {n_assigned:,} / {n_target:,}")
    return clusters


# ==============================================================================
# Step 3: Balanced subsampling (3 cohorts, optionally gated on 10x clusters)
# ==============================================================================
def balanced_subsample(adata, spec: List[Dict[str, Any]],
                        use_official: bool, seed: int):
    np = _import_heavy()["np"]
    banner(3, f"Balanced subsampling (use_official_clusters={use_official})")
    rng = np.random.default_rng(seed)
    chosen: List[str] = []
    for c in spec:
        # Pool: barcodes in the cohort's samples AND in the on-target GEX clusters.
        in_samples = adata.obs["sample_suffix"].isin(c["suffixes"])
        if use_official and "gex_cluster_10x" in adata.obs.columns:
            in_clusters = adata.obs["gex_cluster_10x"].isin(c["gex_clusters"])
            mask = in_samples & in_clusters
            pool_label = "sample_suffix + gex_cluster_10x"
        else:
            mask = in_samples
            pool_label = "sample_suffix (no cluster filter)"
        pool = adata.obs.index[mask]
        n_avail = len(pool)
        n_take  = min(c["n"], n_avail)
        log.info(f"  {c['label']:40s}  pool={pool_label}  "
                 f"available={n_avail:5d}  taking={n_take:5d}")
        if n_take < c["n"]:
            log.warning(f"    Only {n_take} cells match (requested {c['n']}). "
                        f"Downsizing cohort.")
        idx = rng.choice(pool, size=n_take, replace=False)
        chosen.extend(idx)
    chosen = sorted(set(chosen))
    sub = adata[chosen].copy()
    # Map celltypes back onto the subsample.
    celltype_map: Dict[str, str] = {}
    for c in spec:
        in_samples = sub.obs["sample_suffix"].isin(c["suffixes"])
        if use_official and "gex_cluster_10x" in sub.obs.columns:
            in_clusters = sub.obs["gex_cluster_10x"].isin(c["gex_clusters"])
            mask = in_samples & in_clusters
        else:
            mask = in_samples
        celltype_map.update({bc: c["celltype"] for bc in sub.obs.index[mask]})
    sub.obs["CellType"] = sub.obs.index.map(celltype_map)
    n_unmapped = int(sub.obs["CellType"].isna().sum())
    if n_unmapped:
        log.warning(f"  {n_unmapped} cells could not be assigned a CellType "
                    f"(check overlap between cohorts).")
        sub.obs["CellType"] = sub.obs["CellType"].fillna("Unassigned")
    log.info(f"  Final CellType counts:")
    for ct, n in sub.obs["CellType"].value_counts().items():
        log.info(f"    {ct:30s}  {n:4d}")
    return sub


# ==============================================================================
# Step 4: QC report on microglia markers (informational only)
# ==============================================================================
def qc_report(adata, markers: List[str]) -> None:
    banner(4, "QC report on canonical microglia markers")
    np = _import_heavy()["np"]
    sc = _import_heavy()["sc"]
    rna_mask = adata.var["modality"] == "RNA"
    adata_rna = adata[:, rna_mask].copy()
    adata_rna.layers["counts"] = adata_rna.X.copy()
    sc.pp.normalize_total(adata_rna, target_sum=1e4)
    sc.pp.log1p(adata_rna)

    present = [g for g in markers if g in adata_rna.var_names]
    log.info(f"  Markers requested: {markers}")
    log.info(f"  Markers present  : {present}  ({len(present)}/{len(markers)})")
    if len(present) == 0:
        log.warning("  No markers present; skipping the report.")
        return

    X = adata_rna[:, present].X
    if _import_heavy()["sparse"].issparse(X):
        X = X.toarray()
    score = X.mean(axis=1)
    adata.obs["MicrogliaScore"] = score

    log.info(f"  MicrogliaScore summary (per CellType):")
    for ct, vals in adata.obs.groupby("CellType")["MicrogliaScore"]:
        log.info(f"    {ct:30s}  mean={vals.mean():.3f}  "
                 f"median={vals.median():.3f}  min={vals.min():.3f}  max={vals.max():.3f}")

    log.info(f"  Detection rate per marker (fraction of cells with count>0):")
    for g in present:
        if g in adata_rna.var_names:
            col = adata_rna[:, g].X
            if _import_heavy()["sparse"].issparse(col):
                col = col.toarray().ravel()
            frac = (col > 0).mean()
            log.info(f"    {g:10s}  {100*frac:5.1f}%")


# ==============================================================================
# Step 5: Genomic pruning (ATAC peaks to chr7 + chr15 + chr17)
# ==============================================================================
def prune_to_chroms(adata, target_chroms: List[str]) -> None:
    banner(5, f"Pruning ATAC peaks to chromosomes: {target_chroms}")
    atac_mask = adata.var["modality"] == "ATAC"
    chrom_pat = "|".join([f"^{c}[_-]" for c in target_chroms])
    keep = (~atac_mask) | (adata.var_names.to_series().str.match(chrom_pat))
    log.info(f"  ATAC peaks before: {int(atac_mask.sum()):,}")
    log.info(f"  ATAC peaks after : {int((atac_mask & keep).sum()):,}")
    adata._inplace_subset_var(keep)


# ==============================================================================
# Step 6: Barcode synchronization with atac_fragments.tsv.gz
# ==============================================================================
def collect_fragment_barcodes(frag_path: Path) -> set:
    log.info(f"  Scanning {frag_path.name} for unique barcodes...")
    barcodes: set = set()
    n_lines = 0
    t0 = time.time()
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            n_lines += 1
            barcodes.add(line.split("\t", 4)[3])
            if n_lines % 5_000_000 == 0:
                log.info(f"    {n_lines/1e6:6.1f}M fragments scanned "
                         f"({len(barcodes):,} unique barcodes, "
                         f"{time.time()-t0:.1f}s)")
    log.info(f"  Total fragments scanned : {n_lines:,}")
    log.info(f"  Unique barcodes in file : {len(barcodes):,}")
    return barcodes


def synchronize_barcodes(adata, frag_path: Path) -> None:
    banner(6, "Synchronizing barcodes with atac_fragments.tsv.gz")
    frags = collect_fragment_barcodes(frag_path)
    overlap = set(adata.obs_names).intersection(frags)
    n_before = adata.n_obs
    log.info(f"  Subsample barcodes      : {n_before:,}")
    log.info(f"  Overlap with fragments  : {len(overlap):,}  "
             f"({100.0*len(overlap)/max(n_before,1):.1f}%)")
    if len(overlap) < int(0.80 * n_before):
        log.error("Less than 80% barcode overlap with fragments. "
                  "The subsampled cells do not align with the ATAC data. "
                  "This usually means an H5 / fragments mismatch.")
        sys.exit(1)
    adata._inplace_subset_obs(sorted(overlap))
    log.info(f"  Cells after sync: {adata.n_obs:,}")


# ==============================================================================
# Step 7: Filter ATAC fragments to target chromosomes
# ==============================================================================
def filter_fragments_python(frag_in: Path, frag_out: Path,
                            target_chroms: List[str]) -> Tuple[int, int]:
    n_in = n_out = 0
    t0 = time.time()
    chrom_set = set(target_chroms)
    with gzip.open(frag_in, "rt") as fin, \
         gzip.open(frag_out, "wt", compresslevel=6) as fout:
        for line in fin:
            n_in += 1
            chrom = line.split("\t", 1)[0]
            if chrom in chrom_set:
                fout.write(line)
                n_out += 1
            if n_in % 5_000_000 == 0:
                rate = n_in / (time.time() - t0 + 1e-9)
                log.info(f"    {n_in/1e6:6.1f}M scanned, "
                         f"{n_out/1e6:6.1f}M kept "
                         f"({rate/1e6:.1f}M lines/s)")
    return n_in, n_out


def filter_fragments_tabix(frag_in: Path, frag_out: Path,
                           target_chroms: List[str]) -> Tuple[int, int]:
    log.info("  Using system `tabix` for ranged extraction.")
    cmd = ["tabix", "-h", str(frag_in)] + target_chroms
    res = subprocess.run(cmd, check=True, capture_output=True)
    with open(frag_out, "wb") as fh:
        fh.write(res.stdout)
    n_in = 0
    with gzip.open(frag_in, "rt") as f:
        for _ in f:
            n_in += 1
    n_out = 0
    with gzip.open(frag_out, "rt") as f:
        for _ in f:
            n_out += 1
    return n_in, n_out


def index_fragments_tabix(frag_gz: Path) -> None:
    """(Re-)build a .tbi index for the filtered fragments file."""
    if shutil.which("tabix") is None:
        return
    res = subprocess.run(
        ["tabix", "-p", "bed", "-f", str(frag_gz)],
        capture_output=True, text=True,
    )
    if res.returncode == 0:
        log.info(f"  Indexed {frag_gz.name} with tabix.")
    else:
        log.warning(f"  tabix index failed: {res.stderr.strip()}")


def run_fragment_filter(frag_in: Path, frag_out: Path,
                        target_chroms: List[str], no_tabix: bool) -> None:
    banner(7, f"Filtering ATAC fragments to {target_chroms}")
    if not frag_in.is_file() or frag_in.stat().st_size == 0:
        log.error(f"Missing fragments file: {frag_in}")
        sys.exit(1)
    if frag_out.is_file() and frag_out.stat().st_size > 0:
        log.info(f"  Output already exists ({human_size(frag_out.stat().st_size)}). "
                 "Skipping fragment filter. Delete the file to force re-run.")
        return
    have_tabix = shutil.which("tabix") is not None
    if have_tabix and not no_tabix:
        n_in, n_out = filter_fragments_tabix(frag_in, frag_out, target_chroms)
    else:
        if not have_tabix:
            log.info("  `tabix` not on PATH; using pure-Python stream filter.")
        n_in, n_out = filter_fragments_python(frag_in, frag_out, target_chroms)
    log.info(f"  Fragments scanned  : {n_in:,}")
    log.info(f"  Fragments kept     : {n_out:,}  ({100.0*n_out/max(n_in,1):.2f}%)")
    log.info(f"  Output             : {frag_out}  ({human_size(frag_out.stat().st_size)})")
    # (Re-)index for downstream Signac/Cicero usage.
    index_fragments_tabix(frag_out)


# ==============================================================================
# Step 8: Write the .h5ad and (optionally) the Seurat .rds via R shim
# ==============================================================================
def write_h5ad(adata, out_path: Path) -> None:
    banner(8, f"Writing .h5ad (AnnData) -> {out_path.name}")
    adata.write_h5ad(out_path, compression="gzip")
    log.info(f"  OK ({human_size(out_path.stat().st_size)}).")


RDS_SHIM = r"""#!/usr/bin/env Rscript
# Auto-generated by subsample_dataset.py. Converts an .h5ad into a
# SeuratDisk-loaded .rds. Intended to run inside the scMultiome-GRN
# container (which has Seurat, Signac, SeuratDisk installed).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: _write_rds.R <h5ad_in> <rds_out>")
h5ad_in <- args[[1]]
rds_out <- args[[2]]
suppressPackageStartupMessages({
  library(SeuratDisk)
  library(Seurat)
})
tmp_h5seurat <- sub("\\.h5ad$", ".h5seurat", h5ad_in)
Convert(h5ad_in, dest = "h5seurat", overwrite = TRUE, verbose = FALSE)
seu <- LoadH5Seurat(tmp_h5seurat)
if (!"CellType" %in% colnames(seu@meta.data)) {
  warning("CellType column missing after H5Seurat round-trip.")
}
saveRDS(seu, rds_out, compress = "xz")
unlink(tmp_h5seurat)
message(sprintf("Wrote %s (%.2f MB)", rds_out,
                file.info(rds_out)$size / 1024^2))
"""


def write_rds_via_r(h5ad_path: Path, rds_path: Path) -> None:
    banner(8, f"Writing Seurat .rds -> {rds_path.name}")
    if not shutil.which("Rscript"):
        log.warning("  Rscript not on PATH. Skipping .rds generation.")
        log.warning("  Run this script inside the scMultiome-GRN Docker "
                    "container to produce the .rds.")
        return
    shim_path = h5ad_path.parent / "_write_rds.R"
    shim_path.write_text(RDS_SHIM)
    cmd = ["Rscript", str(shim_path), str(h5ad_path), str(rds_path)]
    log.info(f"  Invoking: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        log.error("R shim failed.")
        log.error("STDOUT:\n" + res.stdout)
        log.error("STDERR:\n" + res.stderr)
        log.warning("  Continuing without .rds. The .h5ad is still valid "
                    "and can be converted manually with SeuratDisk::Convert().")
        return
    last = [ln for ln in res.stdout.strip().splitlines() if ln]
    if last:
        log.info("  " + last[-1])
    try:
        shim_path.unlink()
    except OSError:
        pass


# ==============================================================================
# Step 9: Summary JSON
# ==============================================================================
def write_summary(out_dir: Path, summary: Dict[str, Any]) -> None:
    out = out_dir / "subsample_summary.json"
    out.write_text(json.dumps(summary, indent=2, default=str))
    log.info(f"  Summary: {out}")


# ==============================================================================
# Main
# ==============================================================================
def main() -> int:
    args = parse_args()
    t0 = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    _ = _import_heavy()  # fail fast on missing deps

    log.info("")
    log.info("=" * 72)
    log.info(" scMultiome-GRN — Biological test-data subsampler")
    log.info("=" * 72)
    log.info(f"  raw_dir          : {args.raw_dir}")
    log.info(f"  out_dir          : {args.out_dir}")
    log.info(f"  target_chroms    : {args.target_chroms}")
    log.info(f"  seed             : {args.seed}")
    log.info(f"  official_clusters: {args.use_official_10x_clusters}")

    # ---- Validate inputs ----
    log.info("")
    log.info("Validating raw inputs...")
    h5_path       = args.raw_dir / f"{PREFIX}_filtered_feature_bc_matrix.h5"
    frag_path     = args.raw_dir / f"{PREFIX}_atac_fragments.tsv.gz"
    tbi_path      = args.raw_dir / f"{PREFIX}_atac_fragments.tsv.gz.tbi"
    analysis_tar  = args.raw_dir / f"{PREFIX}_analysis.tar.gz"
    require_file(h5_path, "filtered_feature_bc_matrix.h5")
    require_file(frag_path, "atac_fragments.tsv.gz")
    require_file(tbi_path,  "atac_fragments.tsv.gz.tbi")
    if args.use_official_10x_clusters:
        require_file(analysis_tar, "analysis.tar.gz (10x cluster labels)")

    # ---- Resolve the cohort spec ----
    spec = load_cohort_spec(args)

    # ---- Pipeline ----
    adata = load_h5(h5_path)
    annotate_barcodes(adata, spec, args.use_official_10x_clusters)
    if args.use_official_10x_clusters:
        load_10x_clusters(adata, analysis_tar)
    adata_sub = balanced_subsample(adata, spec, args.use_official_10x_clusters, args.seed)
    qc_report(adata_sub, args.microglia_markers)
    prune_to_chroms(adata_sub, args.target_chroms)
    synchronize_barcodes(adata_sub, frag_path)

    frag_out = args.out_dir / f"atac_fragments_{'_'.join(args.target_chroms)}.tsv.gz"
    run_fragment_filter(frag_path, frag_out, args.target_chroms, args.no_tabix)

    h5ad_out = args.out_dir / "microglia_dam_demo.h5ad"
    write_h5ad(adata_sub, h5ad_out)
    if not args.no_rds:
        rds_out = args.out_dir / "microglia_dam_demo.rds"
        write_rds_via_r(h5ad_out, rds_out)

    summary = {
        "n_cells_after_subsample": int(adata_sub.n_obs),
        "celltype_counts": adata_sub.obs["CellType"].value_counts().to_dict(),
        "cohort_spec": spec,
        "n_rna_features":  int((adata_sub.var["modality"] == "RNA").sum()),
        "n_atac_features": int((adata_sub.var["modality"] == "ATAC").sum()),
        "target_chroms":   args.target_chroms,
        "microglia_markers": args.microglia_markers,
        "used_official_clusters": args.use_official_10x_clusters,
        "seed": args.seed,
        "elapsed_seconds": round(time.time() - t0, 2),
    }
    write_summary(args.out_dir, summary)

    log.info("")
    log.info("=" * 72)
    log.info(" DONE")
    log.info("=" * 72)
    log.info(f"  Cells      : {summary['n_cells_after_subsample']}  "
             f"({summary['celltype_counts']})")
    log.info(f"  RNA feats  : {summary['n_rna_features']}")
    log.info(f"  ATAC feats : {summary['n_atac_features']}  "
             f"({' + '.join(args.target_chroms)} only)")
    log.info(f"  Elapsed    : {summary['elapsed_seconds']}s")
    log.info(f"  Outputs in : {args.out_dir}")
    log.info("")
    log.info("Run the Nextflow pipeline with this dataset:")
    log.info("  nextflow run main.nf \\")
    log.info(f"      --input {args.out_dir / 'microglia_dam_demo.rds'} \\")
    log.info("      --genome mm10 \\")
    log.info("      --ko_targets Apoe,Trem2 \\")
    log.info("      -profile docker")
    return 0


if __name__ == "__main__":
    sys.exit(main())
