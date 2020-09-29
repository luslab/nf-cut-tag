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

// Commandline options
/*
    --skip_trim         skip adapter/primer trimming step
    --no_fastqc         skip fastqc


*/


// Define DSL2
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Module global params
-------------------------------------------------------------------------------------------------------------------------------*/
// params {
//     modules {
//         'bowtie2_spike_in'{
//             args             = ""
//             suffix           = ""
//             publish_dir      = "bowtie2_2"
//             publish_results  = "all"
//             unmapped_suffix  = ""
//             output_sam       = false
//         }
//     }
// }

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { luslab_header } from './luslab-nf-modules/tools/luslab_util/main.nf'
include { check_params } from './luslab-nf-modules/tools/luslab_util/main.nf'

include { fastq_metadata } from './luslab-nf-modules/tools/metadata/main.nf'
include { fastqc } from './luslab-nf-modules/tools/fastqc/main.nf'
include { multiqc } from './luslab-nf-modules/tools/multiqc/main.nf'
include { cutadapt } from './luslab-nf-modules/tools/cutadapt/main.nf'
include { bowtie2_align as bt2_genome_align } from './luslab-nf-modules/tools/bowtie2/main.nf'
include { bowtie2_align as bt2_spike_in_align } from './luslab-nf-modules/tools/bowtie2/main.nf'
include { umitools_dedup } from './luslab-nf-modules/tools/umi_tools/main.nf'

// SEACR dev
include { paired_bam_to_bedgraph } from './luslab-nf-modules/workflows/bed_flows/main.nf'

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/
// params.genome_index = ''
// params.spike_in_index = ''

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
    //fastq_metadata.out.view()


    if (!params.no_fastqc){
        // Run fastqc
        fastqc( params.modules['fastqc'], fastq_metadata.out )
        //fastqc.out.view()
    }

    // Run MultiQC?

    // Perform merges here before alignment?

    // Adapter trimming
    cutadapt( params.modules['cutadapt'], fastq_metadata.out )
    //cutadapt.out.fastq.view()

    // Align to genome
    bt2_genome_align( params.modules['bowtie2_align'], cutadapt.out.fastq, params.genome_index )

    // Align to spike-in genome          
    bt2_spike_in_align( params.modules['bowtie2_spike_in'], cutadapt.out.fastq, params.spike_in_index ) 

    // Duplicate removal? 
    //umitools_dedup( params.modules['umi_tools'], bt2_genome_align.out.bam)

    // Spike-in calibration and normalisation

    // Peak-calling
    /// SEACR
    /// prepare BAM files for SEACR input
    paired_bam_to_bedgraph ( bt2_genome_align.out.bam , params.general_genome )

    
    /// Run SEACR



}
