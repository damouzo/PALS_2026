// ==============================================================================
// Module:       prep_r_cicero.nf
// Process:      PREP_R_CICERO
// Description:  R-based Cicero co-accessibility computation over the ATAC
//               assay. Consumes the lightweight Seurat object emitted by
//               PREP_R_SEURAT and emits the all_peaks.csv + cicero_connections
//               .csv tables that the downstream Python INFER_GRN process
//               consumes.
//
//               Runtime: pals-r-bioc (Sequera Wave), R 4.4.3 + Bioconductor
//               3.19 + cicero + monocle3 + EnsDb.Mmusculus.v79. Overridable
//               via --container_r2 <url>.
//
// Inputs:       tuple path(seurat_rds),
//                     path(rna_h5ad),
//                     path(atac_h5ad)        (the seurat_rds is the lightweight
//                                              Seurat object from PREP_R_SEURAT)
//
// Outputs:      tuple path("all_peaks.csv"),
//                     path("cicero_connections.csv"),
//                     path("rna.h5ad"),
//                     path("atac.h5ad")      emit: cicero
//
// Author:       Daniel Mouzo
// License:      MIT
// ==============================================================================

process PREP_R_CICERO {

    tag        "R_cicero_coaccess"
    label      "process_medium"
    publishDir  "${params.outdir}/prep_r_cicero", mode: "copy", overwrite: true

    container  { params.container_r2 }

    input:
    tuple path(seurat_rds),
          path(rna_h5ad),
          path(atac_h5ad)

    output:
    tuple path("all_peaks.csv"),
          path("cicero_connections.csv"),
          path("rna.h5ad"),
          path("atac.h5ad"),
          emit: cicero

    script:
    def args   = task.ext.args   ?: ""
    def seed   = task.ext.seed   ?: params.seed
    def minc   = task.ext.min_count ?: params.min_count
    def maxc   = task.ext.max_count ?: params.max_count
    def mode   = task.ext.mode   ?: params.execution_mode
    """
    set -euo pipefail

    # Run cicero inside the prep_r_cicero publishDir so the CSVs land next to
    # the rna/atac h5ad files copied from the prep_r_seurat step. The h5ad
    # files are already in the workdir thanks to the input staging.
    Rscript ${projectDir}/bin/GRN_cicero.R \\
        --seurat_rds   ${seurat_rds} \\
        --outdir       ./ \\
        --min_count    ${minc} \\
        --max_count    ${maxc} \\
        --seed         ${seed} \\
        --mode         ${mode} \\
        \${args}

    # Sanity: confirm the four expected output files exist.
    for f in all_peaks.csv cicero_connections.csv rna.h5ad atac.h5ad; do
        test -s "\$f" || { echo "ERROR: missing or empty: \$f" >&2; exit 1; }
    done
    """
}
