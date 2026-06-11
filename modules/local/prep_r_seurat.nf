// ==============================================================================
// Module:       prep_r_seurat.nf
// Process:      PREP_R_SEURAT
// Description:  R-based scRNA + scATAC preprocessing with Seurat + Signac +
//               SeuratDisk. Emits the h5ad files and a lightweight Seurat RDS
//               that the downstream PREP_R_CICERO process consumes.
//
//               Runtime: pals-r-seurat (Sequera Wave), R 4.4.3 + Seurat 5.4.0
//               + Signac 1.14.0 + SeuratDisk 0.0.9019 + tidyverse + utility R
//               packages. Overridable via --container_r1 <url>.
//
// Inputs:       path input_rds    Seurat-like object (.rds) with RNA + ATAC
//                                  assays and a CellType meta.data column.
//
// Outputs:      tuple path("rna.h5ad"),
//                     path("atac.h5ad"),
//                     path("seurat_object.rds")    emit: seurat
//
// Author:       Daniel Mouzo
// License:      MIT
// ==============================================================================

process PREP_R_SEURAT {

    tag        "R_seurat_preprocess"
    label      "process_medium"
    publishDir  "${params.outdir}/prep_r_seurat", mode: "copy", overwrite: true

    container  { params.container_r1 }

    input:
    path input_rds

    output:
    tuple path("rna.h5ad"),
          path("atac.h5ad"),
          path("seurat_object.rds"),
          emit: seurat

    script:
    def args   = task.ext.args   ?: ""
    def seed   = task.ext.seed   ?: params.seed
    def minc   = task.ext.min_count ?: params.min_count
    def maxc   = task.ext.max_count ?: params.max_count
    def mode   = task.ext.mode   ?: params.execution_mode
    """
    set -euo pipefail

    Rscript ${projectDir}/bin/GRN_seurat.R \\
        --input_seurat ${input_rds} \\
        --outdir       ./ \\
        --min_count    ${minc} \\
        --max_count    ${maxc} \\
        --seed         ${seed} \\
        --mode         ${mode} \\
        \${args}

    # Sanity: confirm the three expected output files exist.
    for f in rna.h5ad atac.h5ad seurat_object.rds; do
        test -s "\$f" || { echo "ERROR: missing or empty: \$f" >&2; exit 1; }
    done
    """
}
