// ==============================================================================
// Sub-workflow:  multiome_grn.nf
// Purpose:       R -> Python orchestration for the scMultiome GRN pipeline.
//
//                Encapsulates the two-step workflow so it can be included
//                from main.nf (or any other top-level workflow):
//
//                  ch_input (path RDS)
//                       │
//                  [PREPROCESS_ATAC]  (R / Seurat / Signac / Cicero)
//                       │ emits: preproc (all_peaks.csv,
//                       │                     cicero_connections.csv,
//                       │                     rna.h5ad, atac.h5ad)
//                       ▼
//                  [INFER_GRN]        (Python / CellOracle + perturbation)
//
// Inputs:        ch_input  channel emitting the path to the Seurat RDS
//               (other parameters come from the global `params` block)
//
// Emits:         ch_preproc       tuple of 4 emitted paths from PREPROCESS_ATAC
//                ch_grn           tuple of emitted paths from INFER_GRN
//                ch_grn_plots     individual plot paths (for inspection)
//
// Author:        Daniel Mouzo
// License:       MIT
// ==============================================================================

// Sub-workflows need their own process includes (DSL2 scoping).
include { PREPROCESS_ATAC } from '../modules/local/preprocess_atac.nf'
include { INFER_GRN       } from '../modules/local/infer_grn.nf'

workflow GRN_PIPELINE {

    take:
    ch_input                              // path: Seurat RDS file

    main:
    // Step 1: R preprocessing + Cicero co-accessibility
    PREPROCESS_ATAC(ch_input)
    ch_preproc = PREPROCESS_ATAC.out.preproc

    // Step 2: Python CellOracle GRN inference + in silico perturbation
    //
    // The R outputs arrive as a tuple [peaks, connections, rna_h5ad, atac_h5ad].
    // We map the tuple to the 4 path inputs expected by INFER_GRN, and pair
    // them with the seven global parameters the process needs. The order
    // must match the input signature in modules/local/infer_grn.nf
    // (vals first, then paths).
    ch_infer_inputs = ch_preproc
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
    preproc = ch_preproc
    grn     = ch_grn
    plots   = ch_grn_plots
}
