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

    --genome            fasta file of genome
    --bt2_index         bowtie2 index of genome (used instead of genome fasta file)
    --spike_in_genome   spike-in genome, either in fasta/compressed fasta (or bowtie2 index)??


    //TODO add options for overriding module params, especially bt2 aligner
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

include { luslab_header; build_debug_param_summary; check_params } from './luslab-nf-modules/tools/luslab_util/main.nf'
include { fastq_metadata } from './luslab-nf-modules/tools/metadata/main.nf'
include { fastqc } from './luslab-nf-modules/tools/fastqc/main.nf'
include { cutadapt } from './luslab-nf-modules/tools/cutadapt/main.nf'

include { multiqc } from './luslab-nf-modules/tools/multiqc/main.nf'
include { bowtie2_build as bt2_build_exp; bowtie2_build as bt2_build_spike} from './luslab-nf-modules/tools/bowtie2/main.nf'
include { bowtie2_align as bt2_align_exp; bowtie2_align as bt2_align_spike_in } from './luslab-nf-modules/tools/bowtie2/main.nf'

//include { multiqc as multiqc_control} from './luslab-nf-modules/tools/multiqc/main.nf'
//include { cutadapt } from './luslab-nf-modules/tools/cutadapt/main.nf'
//include { bowtie2_align as bt2_genome_align } from './luslab-nf-modules/tools/bowtie2/main.nf'
//include { bowtie2_align as bt2_spike_in_align } from './luslab-nf-modules/tools/bowtie2/main.nf'
//include { umitools_dedup } from './luslab-nf-modules/tools/umi_tools/main.nf'
//include { seacr } from './luslab-nf-modules/tools/seacr/main.nf'

// SEACR dev
//include { paired_bam_to_bedgraph as seacr_data_input} from './luslab-nf-modules/workflows/bed_flows/main.nf'
//include { paired_bam_to_bedgraph as seacr_control_input} from './luslab-nf-modules/workflows/bed_flows/main.nf'
//include { paired_bam_to_bedgraph as seacr_data_control} from './luslab-nf-modules/workflows/bed_flows/main.nf'

/*-----------------------------------------------------------------------------------------------------------------------------
Sub workflows
-------------------------------------------------------------------------------------------------------------------------------*/

//include { qc_align as qc_align_exp; qc_align as qc_align_ctr } from './workflows/qc_align/main.nf'

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/
//params.no_fastqc = ''
// params.input = ''
// params.control = ''
// params.genome_index = ''
// params.spike_in_index = ''

// Module paramaters may need to be different for the workflow concerning the C+T data and the control data
// CUT&Tag data params
//params.cut_tag_params = params.modules

// Control data params 
//params.control_params = params.modules

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/

// Show banner
log.info luslab_header()

// Debug params
if(params.verbose){
    log.info build_debug_param_summary()
}

// Show work summary

// Check inputs
// check_params(['genome']) //will throw up error if these required parameters are not provided 

/*-----------------------------------------------------------------------------------------------------------------------------
Channel setup
-------------------------------------------------------------------------------------------------------------------------------*/
//ch_initial = Channel.from( params.input, params.control )

Channel
    .fromPath("$baseDir/assets/multiqc_config.yml")
    .set { ch_multiqc_config }

if (params.genome) {
    Channel
        .fromPath(params.genome)
        .set { ch_genome }
}

if (params.spike_in_genome) {
    Channel
        .fromPath(params.genome)
        .set { ch_spike_in_genome }
}



// Channel
//     .fromPath( pre_peak_process_data.out.fastqc_path )
//     .mix( pre_peak_process_data.out.cutadapt_path, pre_peak_process_data.out.bt2_path )
//     .collect()
//     .set { ch_data_multiqc }

// Channel
//     .fromPath( pre_peak_process_control.out.fastqc_path )
//     .mix( pre_peak_process_control.out.cutadapt_path, pre_peak_process_control.out.bt2_path )
//     .collect()
//     .set { ch_control_multiqc }

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

workflow {
    /* ---------- Parameter checks to see if bowtie2 indexes need to be built ---------*/

    if (params.genome && !params.bt2_index) {
        bt2_build_exp( params.modules['bowtie2_build'], ch_genome )

        ch_bt2_index = bt2_build_exp.out.bowtieIndex.collect()

    } else {
        Channel
            .fromPath(params.bt2_index)
            .set { ch_bt2_index }
    }

    // TODO Auto-detect spike-in genome 
    // Make bowtie2 index of spike-in genome
    if (params.spike_in_genome) {
        bt2_build_spike( params.modules['bowtie2_build'], ch_spike_in_genome )

        ch_bt2_spike_in = bt2_build_spike.out.bowtieIndex.collect()

    } else {
        Channel
            .fromPath(params.spike_in_genome)
            .set { ch_bt2_spike_in }
    }

    /* ---------- Main Workflow ---------*/

    // Load design file
    fastq_metadata( params.input )

    // Run fastqc
    fastqc( params.modules['fastqc'], fastq_metadata.out.metadata )

    // Adapter trimming
    cutadapt( params.modules['cutadapt'], fastq_metadata.out.metadata )

    // Align to genome
    bt2_align_exp( params.modules['bowtie2_align_exp'], cutadapt.out.fastq, ch_bt2_index )

    bt2_align_exp.out.report | view

    // Align to spike-in genome
    bt2_align_spike_in( params.modules['bowtie2_align_spike_in'], cutadapt.out.fastq, ch_bt2_spike_in )

    bt2_align_spike_in.out.report | view
    
    // Collect reports to produce MultiQC reports
    multiqc( params.modules['multiqc_custom'], ch_multiqc_config, 
        fastqc.out.report
        .mix(cutadapt.out.report)
        .mix(bt2_align_exp.out.report)
        .mix(bt2_align_spike_in.out.report)
        .collect() )


//-------------------------------------------------------------------------------------------------------------------------------*/


    // fastq_metadata.out.metadata | view

    //bowtie2_build( params.modules['bowtie2_build'], ch_genome )

    // Auto-detect index or genome to index
    // if (params.genome && !params.bt2_index) {
    //     bt2_build_exp( params.modules['bowtie2_build'], ch_genome )

    //     ch_bt2_index = bowtie2_build.out.bowtieIndex.collect()

    // } else {
    //     Channel
    //         .fromPath(params.bt2_index)
    //         .set { ch_bt2_index }
    // }

    //     // TODO Auto-detect spike-in genome 
    //     // Make bowtie2 index of spike-in genome
    // if (params.spike_in_genome) {
    //     bt2_build_spike( params.modules['bowtie2_build'], ch_spike_in_genome )

    //     ch_bt2_spike_in = bowtie2_build.out.bowtieIndex.collect()

    // } else {
    //     Channel
    //         .fromPath(params.spike_in_genome)
    //         .set { ch_bt2_spike_in }
    // }

    // ch_bt2_index | view
    // bowtie2_build.out.bowtieIndex | view
    // bowtie2_build.out.report | view

    //qc_align_exp( fastq_metadata.out.metadata, ch_bt2_index, params)
    //qc_align_exp( fastq_metadata.out.metadata, bowtie2_build.out.bowtieIndex.collect(), params)

    //qc_align_exp.out.bam | view
    //qc_align_exp.out.bt2_report | view
    //qc_align_exp.out.fastqc_report | view
    //qc_align_exp.out.cutadapt_report | view

    //qc_align_exp.out.fastqc_report.mix(qc_align_exp.out.cutadapt_report).collect() | view
        //.subscribe {log.info("$it")}

    //TODO Spike-in alignemnet



    // Collect reports to produce MultiQC report

    // multiqc( params.modules['multiqc_custom'], ch_multiqc_config, 
    //     qc_align_exp.out.fastqc_report
    //     .mix(qc_align_exp.out.cutadapt_report)
    //     .mix(qc_align_exp.out.bt2_report)
    //     .collect() )

    // Input data processing
    //pre_peak_process_data( params.input, params.cut_tag_params , params.genome_index )
    //seacr_data_input( pre_peak_process_data.out.bam_file, params.general_genome)

    // Control data processing 
    //pre_peak_process_control( params.control, params.control_params, params.genome_index ) 
    //seacr_control_input( pre_peak_process_control.out.bam_file, params.general_genome)


    // Channels

    //Channel
    //    .value( pre_peak_process_data.out.fastqc_path )
    //    .mix( pre_peak_process_data.out.cutadapt_path, pre_peak_process_data.out.bt2_path )
        //.collect()
    //    .set { ch_data_multiqc }

   // Channel
    //    .value( pre_peak_process_control.out.fastqc_path )
    //    .mix( pre_peak_process_control.out.cutadapt_path, pre_peak_process_control.out.bt2_path )
        //.collect()
    //    .set { ch_control_multiqc }

    // ch_data_multiqc.view()
    // ch_control_multiqc.view()

    // MultiQC - for now implement for each experiment, but see if we can do onComplete too
    // Data MultiQC
   // multiqc_data( params.modules['multiqc'], ch_multiqc_config, ch_data_multiqc )

    // Control MultiQC
   // multiqc_control( params.modules['multiqc'], ch_multiqc_config, ch_control_multiqc )

    // SEACR peak caller
   // seacr( params.modules['seacr'], seacr_data_input.out.bedgraph, seacr_control_input.out.bedgraph )
}




    ////////////////////////////////////// Pre-workflow split up////////////////////////////////////////////////////
    // // Structure input
    // fastq_metadata( params.input )
    // //fastq_metadata.out.view()


    // if (!params.no_fastqc){
    //     // Run fastqc
    //     fastqc( params.modules['fastqc'], fastq_metadata.out )
    //     //fastqc.out.view()
    // }

    // // Run MultiQC?

    // // Perform merges here before alignment?

    // // Adapter trimming
    // cutadapt( params.modules['cutadapt'], fastq_metadata.out )
    // //cutadapt.out.fastq.view()

    // // Align to genome
    // bt2_genome_align( params.modules['bowtie2_align'], cutadapt.out.fastq, params.genome_index )

    // // Align to spike-in genome          
    // //bt2_spike_in_align( params.modules['bowtie2_spike_in'], cutadapt.out.fastq, params.spike_in_index ) 

    // // Duplicate removal? 
    // //umitools_dedup( params.modules['umi_tools'], bt2_genome_align.out.bam)

    // // Spike-in calibration and normalisation

    // // Peak-calling
    // /// SEACR
    // /// prepare BAM files for SEACR input
    // //seacr_data_input ( bt2_genome_align.out.bam , params.general_genome )
    // //seacr_data_contol ( )

    // /// Run SEACR
    // //seacr( params.modules['seacr'], seacr_data_input.out.bedgraph, params.control) 




