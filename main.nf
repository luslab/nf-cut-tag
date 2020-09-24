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
params {
    modules {
        'bowtie2_spike_in'{
            args             = ""
            suffix           = ""
            publish_dir      = "bowtie2"
            publish_results  = "all"
            unmapped_suffix  = ""
            output_sam       = false
        }
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { luslab_header } from './luslab-nf-modules/tools/luslab_util/main.nf'
include { check_params } from './luslab-nf-modules/tools/luslab_util/main.nf'

include { fastq_metadata } from './luslab-nf-modules/tools/metadata/main.nf'
include { fastqc } from './luslab-nf-modules/tools/fastqc/main.nf'
include { cutadapt } from './luslab-nf-modules/tools/cutadapt/main.nf'
include { bowtie2 } from './luslab-nf-modules/tools/bowtie2/main.nf'


/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/
params.genome_index = ''
params.spike_in_index = ''

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
    cutadapt(params.modules['cutadapt'], fastq_metadata.out)
    //cutadapt.out.view()

    // Align to genome
    bowtie2(params.modules['bowtie2'], cutadapt.out, genome_index)

    // Align to spike-in genome
    bowtie2(params.modules['bowtie2_spike_in'], cutadapt.out, spike_in_index)
}
