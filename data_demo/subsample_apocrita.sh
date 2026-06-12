#!/usr/bin/env bash
# ==============================================================================
# data_demo/subsample_apocrita.sh
# ==============================================================================
#
# Apocrita (QMUL HPC) batch script to generate the biological microglia test
# dataset that the scMultiome-GRN pipeline expects at:
#
#   data_demo/processed/microglia_dam_demo.rds
#
# The full Mouse Brain Alzheimer's AppNote multiome matrix is ~33k cells x
# 99k features (sparse, ~1.8 GB in memory) which DOES NOT fit on a typical
# laptop. On Apocrita we run it inside the same Python container the
# pipeline uses (ghcr.io/damouzo/pals-python-grn:1.0.0) and pull the resulting
# .rds back to the laptop for the live demo.
#
# Usage (from the project root, after data_demo/download_data.sh has been
# run so data_demo/raw/ is populated):
#
#   sbatch data_demo/subsample_apocrita.sh
#
# Resources requested:
#   --mem=10GB   (peak RSS of subsample_dataset.py on the full matrix)
#   --cpus=4     (h5py + scanpy + anndata scale linearly up to ~4 threads)
#   --time=2h    (dominated by the 5.8 GB atac_fragments.tsv.gz tabix pass)
#
# The script writes the produced .rds to a host-visible directory
# (defaults to data_demo/processed/) so the next `sbatch` is a no-op
# (idempotent: skips if the .rds already exists and is non-empty).
#
# Apptainer (Singularity) is used instead of Docker because Apocrita nodes
# do not ship with Docker. The image is converted on the fly from the
# Docker image hosted on GHCR.
# ==============================================================================

#SBATCH --job-name=pals-subsample
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -euo pipefail

# ---- Paths -------------------------------------------------------------------
PROJECT_DIR="${PROJECT_DIR:-${PWD}}"
RAW_DIR="${RAW_DIR:-${PROJECT_DIR}/data_demo/raw}"
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/data_demo/processed}"
OUT_RDS="${OUT_RDS:-${OUT_DIR}/microglia_dam_demo.rds}"

# Container reference (matches nextflow.config -> params.container_py).
CONTAINER="docker://ghcr.io/damouzo/pals-python-grn:1.0.0"

mkdir -p "${OUT_DIR}"

# ---- Idempotency check -------------------------------------------------------
if [[ -s "${OUT_RDS}" ]]; then
    echo "[skip] ${OUT_RDS} already exists and is non-empty (size=$(du -h "${OUT_RDS}" | cut -f1))."
    echo "[skip] Delete it manually if you want to regenerate."
    exit 0
fi

# ---- Sanity: raw data must be present ----------------------------------------
if [[ ! -f "${RAW_DIR}/Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_filtered_feature_bc_matrix.h5" ]]; then
    echo "[error] Raw H5 not found at ${RAW_DIR}/."
    echo "[error] Run this first:  bash data_demo/download_data.sh"
    exit 1
fi

# ---- Load Apocrita modules ---------------------------------------------------
module purge
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || {
    echo "[error] Neither 'apptainer' nor 'singularity' module is available on this host."
    exit 1
}

# ---- Pull the container once, then cache it in $HOME ------------------------
APPTAINER_CACHEDIR="${HOME}/.apptainer/cache"
SIF_PATH="${APPTAINER_CACHEDIR}/pals-python-grn_1.0.0.sif"
mkdir -p "${APPTAINER_CACHEDIR}"
if [[ ! -f "${SIF_PATH}" ]]; then
    echo "[info] Pulling and converting ${CONTAINER} -> ${SIF_PATH}"
    apptainer pull --name "${SIF_PATH}" "${CONTAINER}"
fi

# ---- Run the subsample -------------------------------------------------------
echo "[info] Running subsample_dataset.py inside the container"
echo "[info]   project : ${PROJECT_DIR}"
echo "[info]   raw_dir : ${RAW_DIR}"
echo "[info]   out_dir : ${OUT_DIR}"

apptainer exec \
    --bind "${PROJECT_DIR}:/workspace" \
    --pwd   "/workspace" \
    --env   "HOME=/tmp" \
    --env   "MPLCONFIGDIR=/tmp/matplotlib" \
    --env   "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env   "XDG_CACHE_HOME=/tmp/xdg" \
    --env   "PYTHONUNBUFFERED=1" \
    "${SIF_PATH}" \
    python /workspace/data_demo/subsample_dataset.py \
        --raw-dir "/workspace/data_demo/raw" \
        --out-dir "/workspace/data_demo/processed"

# ---- Verify ------------------------------------------------------------------
if [[ ! -s "${OUT_RDS}" ]]; then
    echo "[error] Subsample finished but ${OUT_RDS} is missing or empty."
    echo "[error] Check the SLURM output for tracebacks."
    exit 1
fi

echo "[done] ${OUT_RDS}  (size=$(du -h "${OUT_RDS}" | cut -f1))"
echo "[next] Copy it back to your laptop:"
echo "       scp ${HOSTNAME}:${OUT_RDS} ${OUT_DIR}/"
