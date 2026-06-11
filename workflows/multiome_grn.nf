// ==============================================================================
// Sub-workflow:  multiome_grn.nf
// Purpose:       R -> R -> Python orchestration for the scMultiome GRN pipeline.
//
//                Three-step workflow with one input (Seurat RDS) and a chain of
//                three specialized processes:
//
//                  ch_input (path RDS)
//                       │
//                  [PREP_R_SEURAT]   (R / Seurat + Signac + SeuratDisk, pals-r-seurat)
//                       │ emits: seurat (rna.h5ad, atac.h5ad, seurat_object.rds)
//                       ▼
//                  [PREP_R_CICERO]   (R / cicero + monocle3, pals-r-bioc)
//                       │ emits: cicero (all_peaks.csv, cicero_connections.csv,
//                       │                     rna.h5ad, atac.h5ad)
//                       ▼
//                  [INFER_GRN]       (Python / CellOracle + perturbation, ghcr.io)
//
// Inputs:        ch_input  channel emitting the path to the Seurat RDS
//               (other parameters come from the global `params` block)
//
// Emits:         ch_seurat       tuple of 3 emitted paths from PREP_R_SEURAT
//                ch_cicero       tuple of 4 emitted paths from PREP_R_CICERO
//                ch_grn          tuple of emitted paths from INFER_GRN
//                ch_grn_plots    individual plot paths (for inspection)
//
// Author:        Daniel Mouzo
// License:       MIT
// ==============================================================================

// Sub-workflows need their own process includes (DSL2 scoping).
include { PREP_R_SEURAT } from '../modules/local/prep_r_seurat.nf'
include { PREP_R_CICERO } from '../modules/local/prep_r_cicero.nf'
include { INFER_GRN     } from '../modules/local/infer_grn.nf'

workflow GRN_PIPELINE {

    take:
    ch_input                              // path: Seurat RDS file

    main:
    // Step 1: R Seurat + Signac + SeuratDisk preprocessing
    //   Emits: rna.h5ad, atac.h5ad, seurat_object.rds
    PREP_R_SEURAT(ch_input)
    ch_seurat = PREP_R_SEURAT.out.seurat

    // Step 2: R Cicero co-accessibility (consumes the lightweight Seurat RDS
    //   plus the rna/atac h5ad files emitted by step 1). The h5ad files are
    //   re-emitted so the next step does not need to look them up again.
    ch_cicero_inputs = ch_seurat
        .map { rna_h5ad, atac_h5ad, seurat_rds ->
            tuple(file(seurat_rds), file(rna_h5ad), file(atac_h5ad))
        }
    PREP_R_CICERO(ch_cicero_inputs)
    ch_cicero = PREP_R_CICERO.out.cicero

    // Step 3: Python CellOracle GRN inference + in silico perturbation
    //
    // The R outputs arrive as a tuple [all_peaks.csv, cicero_connections.csv,
    // rna.h5ad, atac.h5ad]. We map the tuple to the 4 path inputs expected by
    // INFER_GRN, and pair them with the seven global parameters the process
    // needs. The order must match the input signature in
    // modules/local/infer_grn.nf (vals first, then paths).
    ch_infer_inputs = ch_cicero
        .map { peaks, connections, rna_h5ad, atac_h5ad ->
            tuple(
                params.percentile,
                params.motif_db,
                params.genome,
                params.clustering,
                params.umap_key,
                params.clust2study,
                params.ko_targets,
                file(peaks), file(connections), file(rna_h5ad), file(atac_h5ad)
            )
        }
    INFER_GRN(ch_infer_inputs)
    ch_grn       = INFER_GRN.out.grn
    ch_grn_plots = INFER_GRN.out.grn.map { _base_grn, _oracle, _links, _scores, plots -> plots }

    emit:
    seurat = ch_seurat
    cicero = ch_cicero
    grn    = ch_grn
    plots  = ch_grn_plots
}
