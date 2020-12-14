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
    
    output:
        path("seq_summary_seaborn.png")
        path("seq_summary_table.csv")

    script:
    """
    python ${meta_script_py} ${meta_table}
    """
}