#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process python_charting {
    publishDir "${params.outdir}/charting",
    mode: "copy", 
    overwrite: true

    container 'quay.io/biocontainers/pybda:0.1.0--pyh5ca1d4c_0'

    input:
        path meta_script_py
        path meta_table
        path(reports)

    output:
        path("alignment_summary.png")
        path("alignment_summary_table.csv")
        path("fragmanet_distribution_violin.png")
        path("fragmanet_distribution_line.png")
        path("fragmanet_distribution_violin.csv")
        path("fragmanet_distribution_line.csv")

    script:
    """
    python ${meta_script_py} ${meta_table} ./
    """
}