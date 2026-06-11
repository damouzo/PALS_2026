// ==============================================================================
// Module:       preprocess_atac.nf
// Process:      PREPROCESS_ATAC
// Description:  R-based ATAC/RNA preprocessing + Cicero co-accessibility.
//
//               Wraps bin/GRN_dataProcess.R and emits four output files:
//                 * all_peaks.csv           (chr-prefixed peaks)
//                 * cicero_connections.csv  (Peak1, Peak2, coaccess)
//                 * rna.h5ad                (SeuratDisk export, RNA assay)
//                 * atac.h5ad               (SeuratDisk export, ATAC assay)
//
// Inputs:       path input_rds    Seurat-like object (.rds)
//
// Outputs:      tuple path("all_peaks.csv"),
//                     path("cicero_connections.csv"),
//                     path("rna.h5ad"),
//                     path("atac.h5ad")        emit: preproc
//
// Author:       Daniel Mouzo
// License:      MIT
// ==============================================================================

process PREPROCESS_ATAC {

    tag        "ATAC_preprocess"
    label      "process_medium"
    publishDir  "${params.outdir}/preprocess_atac", mode: "copy", overwrite: true

    container  { params.container_tag
                  ? "ghcr.io/dmouzo/scmultiome-grn:${params.container_tag}"
                  : "ghcr.io/dmouzo/scmultiome-grn:latest" }

    input:
    path input_rds

    output:
    tuple path("all_peaks.csv"),
          path("cicero_connections.csv"),
          path("rna.h5ad"),
          path("atac.h5ad"),
          emit: preproc

    script:
    def args   = task.ext.args   ?: ""
    def seed   = task.ext.seed   ?: params.seed
    def minc   = task.ext.min_count ?: params.min_count
    def maxc   = task.ext.max_count ?: params.max_count
    """
    set -euo pipefail

    Rscript ${projectDir}/bin/GRN_dataProcess.R \\
        --input_seurat ${input_rds} \\
        --outdir       ./ \\
        --min_count    ${minc} \\
        --max_count    ${maxc} \\
        --seed         ${seed} \\
        ${args}

    # Sanity: confirm the four expected output files exist.
    for f in all_peaks.csv cicero_connections.csv rna.h5ad atac.h5ad; do
        test -s "\$f" || { echo "ERROR: missing or empty: \$f" >&2; exit 1; }
    done
    """
}
