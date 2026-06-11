// ==============================================================================
// Module:       infer_grn.nf
// Process:      INFER_GRN
// Description:  Python-based CellOracle GRN inference + in silico perturbation.
//
//               Wraps bin/GRN_analysis.py. The tag embeds the `--percentile`
//               parameter so changes to it are immediately visible in the
//               task tag and trigger cache invalidation (the -resume demo
//               "wow factor").
//
// Inputs:       tuple path(peaks),
//                     path(connections),
//                     path(rna_h5ad),
//                     path(atac_h5ad),
//                     val(percentile),
//                     val(motif_db),
//                     val(genome),
//                     val(clustering),
//                     val(umap_key),
//                     val(clust2study),
//                     val(ko_targets)
//
// Outputs:      tuple path("base_GRN.parquet"),
//                     path("oracle.celloracle.oracle"),
//                     path("links.celloracle.links"),
//                     path("plots/*.png")     emit: grn
//
// Author:       Daniel Mouzo
// License:      MIT
// ==============================================================================

process INFER_GRN {

    // ⭐ Dynamic tag with the percentile value: any change to --percentile
    // produces a visible tag delta in the Nextflow log AND changes the task
    // hash, which is exactly what makes the -resume cache demo work.
    tag        { "GRN_inference_p${percentile}_${ko_targets.replace(',', '-')}" }
    label      "process_high"
    publishDir  "${params.outdir}/infer_grn", mode: "copy", overwrite: true

    container  { params.container_tag
                  ? "ghcr.io/dmouzo/scmultiome-grn:${params.container_tag}"
                  : "ghcr.io/dmouzo/scmultiome-grn:latest" }

    input:
    tuple val(percentile), val(motif_db), val(genome),
          val(clustering), val(umap_key), val(clust2study), val(ko_targets),
          path(peaks), path(connections), path(rna_h5ad), path(atac_h5ad)

    output:
    tuple path("base_GRN.{parquet,csv}"),
          path("oracle.celloracle.oracle{,.csv}"),
          path("links.celloracle.links{,.csv}"),
          path("network_scores.{xlsx,csv}"),
          path("plots/*.png"),
          emit: grn

    script:
    def args      = task.ext.args ?: ""
    def seed      = task.ext.seed ?: params.seed
    def p_thresh  = task.ext.p_threshold ?: params.p_threshold
    def n_top     = task.ext.n_top_edges ?: params.n_top_edges
    def mode      = task.ext.mode ?: params.execution_mode
    def py        = task.ext.python ?: params.python
    """
    set -euo pipefail

    mkdir -p plots

    ${py} ${projectDir}/bin/GRN_analysis.py \\
        --peaks         ${peaks} \\
        --connections   ${connections} \\
        --rna_h5ad      ${rna_h5ad} \\
        --atac_h5ad     ${atac_h5ad} \\
        --percentile    ${percentile} \\
        --motif_db      ${motif_db} \\
        --genome        ${genome} \\
        --clustering    ${clustering} \\
        --umap_key      ${umap_key} \\
        --clust2study   ${clust2study} \\
        --ko_targets    "${ko_targets}" \\
        --outdir        ./ \\
        --mode          ${mode} \\
        --p_threshold   ${p_thresh} \\
        --n_top_edges   ${n_top} \\
        --seed          ${seed} \\
        ${args}

    # Sanity: at least one vector-field plot must exist.
    n_plots=\$(ls plots/vector_field_KO_*.png 2>/dev/null | wc -l)
    if [ "\$n_plots" -eq 0 ]; then
        echo "ERROR: no vector_field_KO_*.png plots produced" >&2
        exit 1
    fi
    echo "Produced \${n_plots} vector-field plot(s)."
    """
}
