#!/usr/bin/env nextflow

/*
========================================================================================
                         luslab/nf-cut-tag
========================================================================================
Luscombe lab CUT&Tag analysis pipeline.
 #### Homepage / Documentation
 https://github.com/luslab/nf-cut-tag
----------------------------------------------------------------------------------------
*/

// Define DSL2
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Module global params
-------------------------------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { luslab_header } from './luslab-nf-modules/tools/luslab_util/main.nf'
include { check_params } from './luslab-nf-modules/tools/luslab_util/main.nf'

include { fastq_metadata } from './luslab-nf-modules/tools/metadata/main.nf'
include { fastqc } from './luslab-nf-modules/tools/fastqc/main.nf'

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/

// Show banner
log.info luslab_header()

// Show work summary

// Check inputs
// check_params(['genome']) //will throw up error if these required parameters are not provided 

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

// Channel setup

// Run workflow

workflow {

    // Structure input
    fastq_metadata( params.input )
    fastq_metadata.out.view()

    // Run fastqc (change to multiqc, fastqc doesn't unzip)
    fastqc( params.modules['fastqc'], fastq_metadata.out )
    fastqc.out.view()

    // Perform merges here before alignment?

    // Adapter trimming
    
}
