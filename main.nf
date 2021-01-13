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
    --skip_trim             skip adapter/primer trimming step
    --no_fastqc             skip fastqc

    --genome                fasta file of genome
    --bt2_index             bowtie2 index of genome (used instead of genome fasta file)
    --spike_in_genome       spike-in genome, either in fasta/compressed fasta (or bowtie2 index)??
    --genome_blacklist  
    --spike_in_blacklist    optional
    --chrom_sizes           chromasome size data from UCSC e.g (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
    --genome_gtf            GTF for reference genome


    //TODO add options for overriding module params, especially bt2 aligner



    --trim_nextseq INT      INT is q-score cuttoff for nextseq adaptor trimming
*/

// Define DSL2
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Parameter Initialisation
-------------------------------------------------------------------------------------------------------------------------------*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { luslab_header; build_debug_param_summary; check_params } from './luslab-nf-modules/tools/luslab_util/main.nf'
include { fastq_metadata } from './luslab-nf-modules/tools/metadata/main.nf'
include { decompress as decompress_blacklist; decompress as decompress_spike_blacklist; awk as awk_fai; decompress as decompress_genome; } from './luslab-nf-modules/tools/luslab_linux_tools/main.nf'
include { fastqc } from './luslab-nf-modules/tools/fastqc/main.nf'
include { multiqc } from './luslab-nf-modules/tools/multiqc/main.nf'
include { bowtie2_build as bt2_build_exp; bowtie2_build as bt2_build_spike} from './luslab-nf-modules/tools/bowtie2/main.nf'
include { bowtie2_align as bt2_align_exp; bowtie2_align as bt2_align_spike_in } from './luslab-nf-modules/tools/bowtie2/main.nf'
include { meta_report_annotate as meta_annotate_bt2_exp; meta_report_annotate as meta_annotate_bt2_spike; meta_report_annotate as meta_annotate_dt_exp; meta_report_annotate as meta_annotate_dt_spike;} from './luslab-nf-modules/workflows/report_flows/main.nf'
include { deeptools_bam_pe_fragment_size as dt_fragments_exp; deeptools_bam_pe_fragment_size as dt_fragments_spike } from './luslab-nf-modules/tools/deeptools/main.nf'
include { paired_bam_to_bedgraph } from './luslab-nf-modules/workflows/bed_flows/main.nf'
include { seacr } from './luslab-nf-modules/tools/seacr/main.nf'
include { python_charting } from './modules/python_charting/main.nf'

include { samtools_faidx } from './luslab-nf-modules/tools/samtools/main.nf'




//include { multiqc as multiqc_control} from './luslab-nf-modules/tools/multiqc/main.nf'
//include { cutadapt } from './luslab-nf-modules/tools/cutadapt/main.nf'
//include { bowtie2_align as bt2_genome_align } from './luslab-nf-modules/tools/bowtie2/main.nf'
//include { bowtie2_align as bt2_spike_in_align } from './luslab-nf-modules/tools/bowtie2/main.nf'
//include { umitools_dedup } from './luslab-nf-modules/tools/umi_tools/main.nf'


// SEACR dev
//include { paired_bam_to_bedgraph as seacr_data_input} from './luslab-nf-modules/workflows/bed_flows/main.nf'
//include { paired_bam_to_bedgraph as seacr_control_input} from './luslab-nf-modules/workflows/bed_flows/main.nf'
//include { paired_bam_to_bedgraph as seacr_data_control} from './luslab-nf-modules/workflows/bed_flows/main.nf'

// NF-CORE
include { TRIMGALORE } from './nfcore-nf-modules/software/trimgalore/main' addParams( options: trimgalore_options )

/*-----------------------------------------------------------------------------------------------------------------------------
Sub workflows
-------------------------------------------------------------------------------------------------------------------------------*/

//include { qc_align as qc_align_exp; qc_align as qc_align_ctr } from './workflows/qc_align/main.nf'

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/

// Show banner
log.info luslab_header()

// Debug params
if(params.verbose){
    log.info build_debug_param_summary()
}

// Function for checking file extensions
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Show work summary

// Check inputs
// check_params(['genome']) //will throw up error if these required parameters are not provided 

/*-----------------------------------------------------------------------------------------------------------------------------
Channel setup
-------------------------------------------------------------------------------------------------------------------------------*/

Channel
    .fromPath("$baseDir/assets/multiqc_config.yml")
    .set { ch_multiqc_config }

Channel
    .fromPath(params.chrom_sizes)
    .set { ch_chrom_sizes }

Channel
    .fromPath(params.genome_gtf)
    .set { ch_genome_gtf }

if (params.genome) {
    Channel
        .fromPath(params.genome)
        .set { ch_genome }
}

// Deal with compressed vs decompressed genome fasta file
if (hasExtension(params.genome, 'gz')) {
    meta_genome =[
        [[:], params.genome]
    ]
    Channel
        .from(meta_genome)
        .set { ch_genome_decompress }
} else {
    Channel
        .from(params.genome)
        .set { ch_decompressed_genome }
}

if (params.spike_in_genome) {
    Channel
        .fromPath(params.spike_in_genome)
        .set { ch_spike_in_genome }
}

// Deal with compressed vs decompressed genome blacklist
if (hasExtension(params.genome_blacklist, 'gz')) {
    genome_blacklist = [
        [[:], params.genome_blacklist]
    ]
    Channel
        .from(genome_blacklist)
        .set { ch_genome_blacklist_decompress }
} else {
    Channel
        .from(params.genome_blacklist)
        .set { ch_decompressed_genome_blacklist }
}

// Deal with compressed vs decompressed spike-in blacklist
if (hasExtension(params.spike_in_blacklist, 'gz')) {
    spike_in_blacklist = [
        [[:], params.spike_in_blacklist]
    ]
    Channel
        .from(spike_in_blacklist)
        .set { ch_spike_blacklist_decompress }
} else {
    Channel
        .from(params.spike_in_blacklist)
        .set { ch_decompressed_spike_blacklist }
}

Channel
    .value("$baseDir/assets/awk_scripts/bt2_report_to_csv.awk")
    .set { ch_bt2_awk }

Channel
    .value("$baseDir/assets/awk_scripts/bt2_spike_report_to_csv.awk")
    .set { ch_bt2_spike_awk }

Channel
    .value(params.normalisation_c)
    .set{ ch_normalisation_c }

Channel
    .fromPath("$baseDir/bin/cut_tag_figs.py")
    .set{ ch_charting_script }

Channel
    .value("$baseDir/assets/awk_scripts/dt_report_annotate.awk")
    .set{ ch_dt_awk }

Channel
    .value("$baseDir/assets/awk_scripts/dt_report_annotate.awk")
    .set{ ch_dt_spike_awk }

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
        bt2_build_exp( modules['bowtie2_build'], ch_genome )

        ch_bt2_index = bt2_build_exp.out.bowtieIndex.collect()
    } else {
        Channel
            .fromPath(params.bt2_index)
            .set { ch_bt2_index }
    }

    // TODO Auto-detect spike-in genome 
    // Make bowtie2 index of spike-in genome
    if (params.spike_in_genome) {
        bt2_build_spike( modules['bowtie2_build'], ch_spike_in_genome )

        ch_bt2_spike_in = bt2_build_spike.out.bowtieIndex.collect()

    } else {
        Channel
            .fromPath(params.spike_in_genome)
            .set { ch_bt2_spike_in }
    }

    /* ---------- Decompress files if necessary ---------*/    

    if (hasExtension(params.genome_blacklist, 'gz')) {
        decompress_blacklist ( ch_genome_blacklist_decompress )
        ch_decompressed_genome_blacklist = decompress_blacklist.out.file_no_meta
    }

    if (hasExtension(params.spike_in_blacklist, 'gz')) {
        decompress_spike_blacklist ( ch_spike_blacklist_decompress )
        ch_decompressed_spike_blacklist = decompress_spike_blacklist.out.file_no_meta
    }

    if (hasExtension(params.spike_in_blacklist, 'gz')) {
        decompress_spike_blacklist ( ch_spike_blacklist_decompress )
        ch_decompressed_spike_blacklist = decompress_spike_blacklist.out.file_no_meta
    }

    if (hasExtension(params.genome, 'gz')) {
        decompress_genome( ch_genome_decompress )
        ch_decompressed_genome = decompress_genome.out.file_no_meta
    }

    /* ---------- Main Workflow ---------*/

    // ***** LOAD DESIGN ***** //
    fastq_metadata( params.input )
    //fastq_metadata.out.metadata | view

    // ***** RAW FASTQC ***** //
    fastqc( modules['fastqc'], fastq_metadata.out.metadata )

    // ***** TRIM GALORE (Adaptor trimming and post-trim fastqc) ***** //
    TRIMGALORE ( fastq_metadata.out.metadata )
    //TRIMGALORE.out.reads | view

    // ***** BOWTIE2 - ALIGN TO REFERENCE GENOME + SORT AND INDEX ***** //
    bt2_align_exp( modules['bowtie2_align_exp'], TRIMGALORE.out.reads, ch_bt2_index )
    //bt2_align_exp.out.bam | view

    // ***** ANNOTATE METADATA WITH REFERENCE GENOME ALIGNMENT BT2 STATS ***** //
    meta_annotate_bt2_exp( bt2_align_exp.out.report_meta, bt2_align_exp.out.bam, ch_bt2_awk, modules )
    //meta_annotate_bt2_exp.out.annotated_input | view

    // ***** BOWTIE2 - ALIGN TO SPIKE-IN GENOME + SORT AND INDEX ***** //
    bt2_align_spike_in( modules['bowtie2_align_spike_in'], TRIMGALORE.out.reads, ch_bt2_spike_in )
    //bt2_align_spike_in.out.bam | view

    // ***** ANNOTATE METADATA WITH SPIKE-IN GENOME ALIGNMENT BT2 STATS ***** //
    meta_annotate_bt2_spike( bt2_align_spike_in.out.report_meta, bt2_align_spike_in.out.bam, ch_bt2_spike_awk, params.modules)
    //meta_annotate_bt2_spike.out.annotated_input | view

    // ***** CALCULATE REFERENCE GENOME ALIGNMENT FRAGMENT STATS ***** //
    dt_fragments_exp( modules['deeptools_bam_pe_fragment_size'], meta_annotate_bt2_exp.out.annotated_input, ch_decompressed_genome_blacklist.collect() )
    //dt_fragments_exp.out.fragment_size_summary | view
    
    // ***** ANNOTATE METADATA WITH REFERENCE GENOME ALIGNMENT FRAGMENT STATS ***** //
    meta_annotate_dt_exp( dt_fragments_exp.out.fragment_stats_meta, meta_annotate_bt2_exp.out.annotated_input, ch_dt_awk, params.modules )
    //meta_annotate_dt_exp.out.annotated_input | view

    // ***** CALCULATE SPIKE-IN GENOME ALIGNMENT FRAGMENT STATS ***** //
    dt_fragments_spike( modules['deeptools_bam_pe_fragment_size'], meta_annotate_bt2_spike.out.annotated_input, ch_decompressed_spike_blacklist.collect() )
    //dt_fragments_spike.out.fragment_size_summary | view

    // ***** ANNOTATE METADATA WITH REFERENCE GENOME ALIGNMENT FRAGMENT STATS ***** //
    meta_annotate_dt_spike( dt_fragments_spike.out.fragment_stats_meta, meta_annotate_bt2_spike.out.annotated_input, ch_dt_spike_awk, params.modules )
    //meta_annotate_dt_spike.out.annotated_input | view

    // ***** CHANNEL NAME CLEAN-UP ***** //
    final_meta_exp = meta_annotate_dt_exp.out.annotated_input
    final_meta_spike = meta_annotate_dt_spike.out.annotated_input

    // ***** CALCULATE SCALE FACTOR FOR SPIKE-IN NORMALISATION ***** //
    if (params.spike_in_genome){
        final_meta_spike
            .combine ( ch_normalisation_c )
            .map { row -> [ row[0].sample_id, row[3] / (row[0].find{ it.key == "bt2_spike_total_aligned" }?.value.toInteger()) ] }
            .set { ch_scale_factor }
    } else { // this else doesn't make sense because there would be no spike_in_meta_out from alignment if now spike-in genome is provided
        //spike_in_meta_annotate.out.annotated_input
        final_meta_spike
        .map { row -> [ row[0].sample_id, 1] }
        .set { ch_scale_factor }
    }
    //ch_scale_factor | view

    // ***** MERGE SCALE FACTOR INTO DATA CHANNELS ***** //
    bt2_align_exp.out.bam
        .map { row -> [row[0].sample_id, row ].flatten()}
        .join ( ch_scale_factor )
        .map { row -> row[1..(row.size() - 1)] }
        .multiMap { it ->
            bt2_bam_tuple: it [0..-2]
            scale_factor: it[-1]
        }
        .set { ch_align_scale }
    //ch_align_scale.bt2_bam_tuple.collect() | view
    //ch_align_scale.scale_factor.collect() | view

    // ***** CONVERT DATA BAM's TO BEDGRAPH AND NORMALISE ***** //
    paired_bam_to_bedgraph( ch_align_scale.bt2_bam_tuple, ch_align_scale.scale_factor )
    //paired_bam_to_bedgraph.out.bedgraph | view

    // ***** SPLIT EXPERIMENT AND CTRL DATA ***** //
    paired_bam_to_bedgraph.out.bedgraph
        //.map { row -> [row[0].control, row ].flatten()}
        .branch { it ->
            ch_exp: it[0].control == 'no'
            ch_control: it[0].control == 'yes'
        }
        .set { ch_split }
    //ch_split.ch_exp.collect() | view
    //ch_split.ch_control.collect() | view

    // ***** EXTRACT GROUP FROM META ***** //
    ch_split.ch_control
        .map { row -> [row[0].group, row ].flatten() }
        .set { ch_control_group }
    //ch_control_group | view

    ch_split.ch_exp
        .map { row -> [row[0].group, row ].flatten() }
        .set { ch_exp_group }
    //ch_exp_group | view

    // ***** ORDER SO REPLICATES MATCH UP ***** //
    ch_control_group
        .cross ( ch_exp_group )
        .multiMap { it ->
            ch_exp_bedgraph: it[1][1..-1]
            ch_control_bedgraph: it[0][1..-1]
        }
        .set { ch_exp_ctrl_split }
    //ch_exp_ctrl_split.ch_exp_bedgraph | view
    //ch_exp_ctrl_split.ch_control_bedgraph | view

    // ***** CALL PEAKS ***** //
    seacr( modules['seacr'], ch_exp_ctrl_split.ch_exp_bedgraph, ch_exp_ctrl_split.ch_control_bedgraph )
    //seacr.out.bed | view
    
    // ***** COLLECT MULTIQC AND RUN ***** //
    multiqc( params.modules['multiqc_custom'], ch_multiqc_config, 
        fastqc.out.report
        .mix(TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]))
        .mix(TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]))
        .mix(bt2_align_exp.out.report)
        .mix(bt2_align_spike_in.out.report)
        .collect() )

    // ***** EXTRACT SAMPLE ID FROM META ***** //
    final_meta_exp
        .map { row -> [row[0].sample_id, row[0]].flatten() }
        .set { ch_exp_meta_sample_id }
    //ch_exp_meta_sample_id | view

    final_meta_spike
        .map { row -> [row[0].sample_id, row[0]].flatten() }
        .set { ch_spike_in_meta_sample_id }
    //ch_spike_in_meta_sample_id | view

    // ***** JOIN DATA ON SAMPLE ID ***** //
    ch_exp_meta_sample_id
        .join ( ch_spike_in_meta_sample_id )
        .map { row -> row[1] << row[2] }// [bt2_spike_align1] } //, row[2].find{ it.key == "bt2_spike_align_gt1" }, row[2].find{ it.key == "bt2_spike_non_aligned" }, row[2].find{ it.key == "bt2_spike_total_aligned" } ] }
        .collect()
        .set { ch_meta_all }
    //ch_meta_all | view

    // ***** WRITE META DATA TO TABLE ***** //
    meta_file ( ch_meta_all )
    //meta_file.out.meta_table | view

    // ***** CREATE PLOTS ***** //
    python_charting ( ch_charting_script, meta_file.out.meta_table, dt_fragments_exp.out.fragment_no_meta.collect() )

    // ***** CONVERT BEDGRAPH TO BIGWIG ***** //
    ucsc_bedgraphtobigwig( paired_bam_to_bedgraph.out.bedgraph, ch_chrom_sizes.collect() )

    // ***** CREATE IGV SESSION ***** //
    igv( ch_decompressed_genome.collect(), ch_genome_gtf, seacr.out.bed.collect{it[1]}.ifEmpty([]), ucsc_bedgraphtobigwig.out.bigwig.collect() )
}

process meta_file {
    publishDir "${params.outdir}/meta",
    mode: "copy",
    overwrite: true
    
    container 'ubuntu:16.04'

    input:
        val(all_meta)

    output:
        path "meta_table.csv", emit: meta_table

    script:
    arr_str = all_meta[0].keySet().join(",") + ","

    for ( int i = 0;i<all_meta.size();i++ ) {
        sample_str = all_meta[i].values().join(",")
        arr_str =  arr_str + "\n" + sample_str
    }

    """
    echo "${arr_str}" > meta_table.csv

    """
}

process ucsc_bedgraphtobigwig {
    publishDir "${params.outdir}/bigwig",
    mode: "copy",
    overwrite: true

    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1'

    input:
      tuple val(meta), path(bedgraph)
      path chrome_sizes

    output:
      path "*.bw", emit: bigwig

    script:
    """
    bedGraphToBigWig ${bedgraph} ${chrome_sizes} ${bedgraph.simpleName}.bw
    """
}

process igv {
    publishDir "${params.outdir}/igv",
    mode: "copy",
    overwrite: true

    container 'python:3.9.1-buster'

    input:
      val(genome)
      path gtf
      path beds
      path bw

    output:
      path('*.{txt,xml,bed,bw,gz}', includeInputs:true)

    script:
    """
    find -L * -iname "*.gz" -exec echo -e {}"\\t0,0,178" \\; > gz.igv.txt
    find -L * -iname "*.bed" -exec echo -e {}"\\t0,0,178" \\; > bed.igv.txt
    find -L * -iname "*.bw" -exec echo -e {}"\\t0,0,178" \\; > bw.igv.txt
    cat *.txt > igv_files.txt
    ${baseDir}/bin/igv_files_to_session.py igv_session.xml igv_files.txt ${genome[0]} --path_prefix './'
    """
}

/*--------------------------archive channel manipulations-----------------------------*/
//     // Get scale factor for normalisation
//     if (params.spike_in_genome){
//         spike_in_meta_annotate.out.annotated_input
//             .combine ( ch_normalisation_c )
//             .map { row -> [ row[0].sample_id, row[3] / (row[0].find{ it.key == "bt2_spike_total_aligned" }?.value.toInteger()) ] }
//             .set { ch_scale_factor }
//        // ch_scale_factor | view
//        log.info "got here"
//     } else { // this else doesn't make sense because there would be no spike_in_meta_out from alignment if now spike-in genome is provided
//         spike_in_meta_annotate.out.annotated_input
//         .map { row -> [ row[0].sample_id, 1] }
//         .set { ch_scale_factor }
//     }
//     ch_scale_factor | view
//    // bt2_align_exp.out.bam | view

//     // Align scale factor and sample to parse to paired_bam_to_bedgraph
//     bt2_align_exp.out.bam
//         .map { row -> [row[0].sample_id, row ].flatten()}
//         .join ( ch_scale_factor )
//         .map { row -> row[1..(row.size() - 1)] }
//         //.set { ch_bt2_align_scale }
//         .multiMap { it ->
//             bt2_bam_tuple: it [0..-2]
//             scale_factor: it[-1]
//         }
//         .set { ch_align_scale }
//     //ch_align_scale.bt2_bam_tuple | view
//     // ch_align_scale.scale_factor | view
//     // ch_bt2_align_scale | view

//     // Genome size index not currently needed, commenting out this section
//     // // Produce genome size index
//     // if (hasExtension(params.genome, 'gz')) {
//     //     decompress( ch_genome_decompress )
//     //     ch_decompressed_genome = decompress.out.file_no_meta
//     // }
//     // //decompress.out.file_no_meta | view
//     // samtools_faidx( params.modules['samtools_faidx'], ch_decompressed_genome )
//     // //samtools_faidx.out.fai | view 
//     // awk_fai( params.modules['awk_fai'], samtools_faidx.out.fasta )
//     // // awk_fai.out.file_no_meta | view

//     // Convert bam files to bedgraphs (does not need to be performed on spike-in alignment?)
//     //paired_bam_to_bedgraph( ch_align_scale.bt2_bam_tuple, awk_fai.out.file_no_meta.collect(), ch_align_scale.scale_factor )
//     paired_bam_to_bedgraph( ch_align_scale.bt2_bam_tuple, ch_align_scale.scale_factor )

//     // Split experiment and control
//     paired_bam_to_bedgraph.out.bedgraph
//         //.map { row -> [row[0].control, row ].flatten()}
//         .branch { it ->
//             ch_exp: it[0].control == 'no'
//             ch_control: it[0].control == 'yes'
//         }
//         .set { ch_split }
//     //ch_split.ch_exp | view
//     //ch_split.ch_control | view

//     ch_split.ch_control
//         .map { row -> [row[0].group, row ].flatten() }
//         .set { ch_control_group }
//     // ch_control_group | view

//     ch_split.ch_exp
//         .map { row -> [row[0].group, row ].flatten() }
//         .set { ch_exp_group }
//     // ch_exp_group | view

//     ch_control_group
//         .cross ( ch_exp_group )
//         .multiMap { it ->
//             ch_exp_bedgraph: it[1][1..-1]
//             ch_control_bedgraph: it[0][1..-1]
//         }
//         .set { ch_exp_ctrl_split }
//     // ch_exp_ctrl_split | view
//     // ch_exp_ctrl_split.ch_exp_bedgraph | view
//     // ch_exp_ctrl_split.ch_control_bedgraph | view

//      Then seacr and following stuff
/*--------------------------------------------------------------*/

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




