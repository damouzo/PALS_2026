# ==============================================================================
# data_demo/download_data.sh
#
# Production-grade, idempotent downloader for the 10x Genomics scMultiome
# Mouse Brain Alzheimer's AppNote dataset.
#
# Reference landing page:
#   https://www.10xgenomics.com/datasets/multiomic-integration-neuroscience-application-note-single-cell-multiome-rna-atac-alzheimers-disease-mouse-model-brain-coronal-sections-from-one-hemisphere-over-a-time-course-1-standard
#
# This script fetches the four assets required by data_demo/subsample_dataset.py:
#
#   1. filtered_feature_bc_matrix.h5        (~225 MB) — aggregated RNA + ATAC
#   2. atac_fragments.tsv.gz                (~6.26 GB) — chromatin cuts
#   3. atac_fragments.tsv.gz.tbi            (~2 MB)    — Tabix index
#   4. Multiome_..._web_summary.html        (~10.7 MB) — QC web summary
#
# Idempotency: files that already exist and are non-empty are skipped
# (so re-running the script is a no-op). The largest file uses `curl -C -`
# to support HTTP resume if a previous run was interrupted.
#
# Usage:
#   bash data_demo/download_data.sh
#   bash data_demo/download_data.sh --force        # re-download everything
#   bash data_demo/download_data.sh --help
#
# Requirements: curl OR wget, ~7 GB free disk space.
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------
readonly BASE_URL="https://cf.10xgenomics.com/samples/cell-arc/2.0.1/Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote"
readonly PREFIX="Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly TARGET_DIR="${SCRIPT_DIR}/raw"

# Asset catalog:  filename | expected_size_bytes (used only for logging).
# Sizes are from the 10x Genomics landing page (October 2024 snapshot).
declare -a ASSETS=(
    "${PREFIX}_filtered_feature_bc_matrix.h5|225887713"
    "${PREFIX}_atac_fragments.tsv.gz|6262747100"
    "${PREFIX}_atac_fragments.tsv.gz.tbi|1951093"
    "${PREFIX}_web_summary.html|10751543"
)

# ------------------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------------------
FORCE=0
for arg in "$@"; do
    case "${arg}" in
        --force) FORCE=1 ;;
        --help|-h)
            sed -n '2,30p' "${BASH_SOURCE[0]}"
            exit 0
            ;;
        *)
            printf 'ERROR: unknown argument: %s\n' "${arg}" >&2
            exit 64
            ;;
    esac
done

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
have_curl=0
have_wget=0
command -v curl >/dev/null 2>&1 && have_curl=1
command -v wget >/dev/null 2>&1 && have_wget=1
if [ "${have_curl}" -eq 0 ] && [ "${have_wget}" -eq 0 ]; then
    printf 'ERROR: neither curl nor wget is installed.\n' >&2
    exit 69
fi

log_step() {
    # log_step <current> <total> <message>
    printf '\n[%s/%s] %s\n' "${1}" "${2}" "${3}"
}

log_info() {
    printf '       %s\n' "${1}"
}

log_warn() {
    printf '       [WARN] %s\n' "${1}" >&2
}

log_err() {
    printf '       [ERROR] %s\n' "${1}" >&2
}

human_size() {
    # human_size <bytes>
    local bytes="${1}"
    if   [ "${bytes}" -ge 1073741824 ]; then printf '%.2f GB' "$(echo "${bytes} / 1073741824" | bc -l)"
    elif [ "${bytes}" -ge 1048576 ];    then printf '%.2f MB' "$(echo "${bytes} / 1048576"    | bc -l)"
    elif [ "${bytes}" -ge 1024 ];       then printf '%.2f KB' "$(echo "${bytes} / 1024"       | bc -l)"
    else printf '%d B' "${bytes}"
    fi
}

fetch() {
    # fetch <url> <out_path> <is_resumable>  (sets global FETCHED_BYTES)
    local url="${1}"
    local out="${2}"
    local resumable="${3}"
    FETCHED_BYTES=0
    if [ "${have_curl}" -eq 1 ]; then
        if [ "${resumable}" -eq 1 ]; then
            # -C - enables HTTP resume; -f fails on HTTP errors; -L follows
            # redirects; -sS silent except errors; --retry 3 with backoff.
            curl -fL -C - -sS --retry 3 --retry-delay 5 -o "${out}" "${url}"
        else
            curl -fL -sS --retry 3 --retry-delay 5 -o "${out}" "${url}"
        fi
    else
        if [ "${resumable}" -eq 1 ]; then
            wget -c -q --tries=3 --retry-connrefused --waitretry=5 -O "${out}" "${url}"
        else
            wget -q --tries=3 --retry-connrefused --waitretry=5 -O "${out}" "${url}"
        fi
    fi
    if [ -f "${out}" ]; then
        FETCHED_BYTES=$(stat -c%s "${out}" 2>/dev/null || stat -f%z "${out}")
    fi
}

# ------------------------------------------------------------------------------
# Pre-flight
# ------------------------------------------------------------------------------
printf '\n===========================================================\n'
printf ' scMultiome-GRN — 10x Genomics dataset downloader\n'
printf '===========================================================\n'
printf ' Target dir : %s\n' "${TARGET_DIR}"
printf ' Base URL   : %s\n' "${BASE_URL}"
printf ' HTTP client: %s\n' "$([ "${have_curl}" -eq 1 ] && echo 'curl' || echo 'wget')"
printf ' Force mode : %s\n' "$([ "${FORCE}" -eq 1 ] && echo 'YES (re-download all)' || echo 'NO (skip existing)')"
printf '%s' '-----------------------------------------------------------'
printf '\n'

mkdir -p "${TARGET_DIR}"

START_TS=$(date +%s)
TOTAL=${#ASSETS[@]}
CURRENT=0
TOTAL_BYTES=0
SKIPPED=0
DOWNLOADED=0

for entry in "${ASSETS[@]}"; do
    CURRENT=$((CURRENT + 1))
    fname="${entry%%|*}"
    expected_size="${entry##*|}"
    url="${BASE_URL}/${fname}"
    out_path="${TARGET_DIR}/${fname}"

    log_step "${CURRENT}" "${TOTAL}" "Asset: ${fname}  (expected $(human_size "${expected_size}"))"

    # Skip if already present and non-empty (unless --force).
    if [ -s "${out_path}" ] && [ "${FORCE}" -eq 0 ]; then
        local_size=$(stat -c%s "${out_path}" 2>/dev/null || stat -f%z "${out_path}")
        log_info "Already present (${local_size} bytes). Skipping. Use --force to re-download."
        SKIPPED=$((SKIPPED + 1))
        TOTAL_BYTES=$((TOTAL_BYTES + local_size))
        continue
    fi

    # Resume supported for the large .tsv.gz; not for the rest (no real benefit).
    resumable=0
    case "${fname}" in
        *.tsv.gz) resumable=1 ;;
    esac

    log_info "Fetching: ${url}"
    if fetch "${url}" "${out_path}" "${resumable}"; then
        if [ ! -s "${out_path}" ]; then
            log_err "Download produced an empty file: ${out_path}"
            exit 70
        fi
        log_info "OK: $(human_size "${FETCHED_BYTES}") written to ${out_path}"
        DOWNLOADED=$((DOWNLOADED + 1))
        TOTAL_BYTES=$((TOTAL_BYTES + FETCHED_BYTES))
    else
        log_err "Download failed: ${url}"
        exit 70
    fi
done

# ------------------------------------------------------------------------------
# Post-flight validation
# ------------------------------------------------------------------------------
printf '\n%s' '-----------------------------------------------------------'
printf '\n'
printf ' Post-download validation\n'
printf '%s' '-----------------------------------------------------------'
printf '\n'
all_ok=1
for entry in "${ASSETS[@]}"; do
    fname="${entry%%|*}"
    out_path="${TARGET_DIR}/${fname}"
    if [ ! -s "${out_path}" ]; then
        log_err "Missing or empty: ${fname}"
        all_ok=0
    else
        sz=$(stat -c%s "${out_path}" 2>/dev/null || stat -f%z "${out_path}")
        log_info "  [OK]   ${fname}  ($(human_size "${sz}"))"
    fi
done

if [ "${all_ok}" -ne 1 ]; then
    log_err "One or more assets are missing or empty. Re-run the script (it resumes)."
    exit 70
fi

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))

printf '\n===========================================================\n'
printf ' Download summary\n'
printf '===========================================================\n'
printf ' Downloaded : %d file(s)\n' "${DOWNLOADED}"
printf ' Skipped    : %d file(s)\n' "${SKIPPED}"
printf ' Total size : %s\n' "$(human_size "${TOTAL_BYTES}")"
printf ' Elapsed    : %d second(s)\n' "${ELAPSED}"
printf ' Target     : %s\n' "${TARGET_DIR}"
printf '%s' '-----------------------------------------------------------'
printf '\n'
printf ' Next step  : python data_demo/subsample_dataset.py\n'
printf '===========================================================\n\n'
