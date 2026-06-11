#!/usr/bin/env nextflow

// ==============================================================================
// scMultiome-GRN — Pipeline entry point
// ==============================================================================
//
// Native polyglot Nextflow DSL2 pipeline that integrates scRNA-seq + scATAC-seq
// data, infers a base Gene Regulatory Network (CellOracle) and runs in silico
// transcription-factor perturbation simulations.
//
//   R stage  :  Seurat / Signac / Cicero     (PREPROCESS_ATAC)
//   Python   :  CellOracle / GimmeMotifs    (INFER_GRN)
//
// Run with the bundled test dataset and Docker:
//
//   nextflow run main.nf -profile test,docker
//
// Re-run with a different network-pruning percentile (the -resume demo):
//
//   nextflow run main.nf -profile test,docker --percentile 95 -resume
//
// ==============================================================================
//
// All default parameters are declared in `nextflow.config` (see conf/ for
// resource profiles and modules.config for per-process overrides).
//
// ==============================================================================

nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// Module / sub-workflow includes
// -----------------------------------------------------------------------------
include { GRN_PIPELINE } from './workflows/multiome_grn.nf'

// -----------------------------------------------------------------------------
// Main workflow
// -----------------------------------------------------------------------------
workflow {

    log.info """
    ===========================================================
     scMultiome-GRN pipeline
    ===========================================================
      input          : ${params.input}
      outdir         : ${params.outdir}
      percentile     : ${params.percentile}
      ko_targets     : ${params.ko_targets}
      genome         : ${params.genome}
      mode           : ${params.execution_mode}
    ===========================================================
    """.stripIndent()

    // Build the input channel from a single file path.
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .map { f -> tuple(f) }

    // Run the polyglot pipeline via the reusable sub-workflow.
    GRN_PIPELINE(ch_input)

    // Final reporting on what was emitted.
    GRN_PIPELINE.out.plots.view { plot_path ->
        log.info "Plot produced: ${plot_path}"
    }
}
