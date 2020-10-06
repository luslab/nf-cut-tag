#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Module inclusions

include { fastq_metadata } from '../luslab-nf-modules/tools/metadata/main.nf'
include { fastqc } from '../luslab-nf-modules/tools/fastqc/main.nf'
include { cutadapt } from '../luslab-nf-modules/tools/cutadapt/main.nf'
include { bowtie2_align } from '../luslab-nf-modules/tools/bowtie2/main.nf'

workflow pre_peak_process {
    take: reads
    take: flow_params
    take: genome_ind
    main:

    // Structure input
    fastq_metadata( reads )

    // Run fastqc
    fastqc( flow_params['fastqc'], fastq_metadata.out )

    // Adapter trimming
    cutadapt( flow_params['cutadapt'], fastq_metadata.out )

    // Align to genome
    bowtie2_align( flow_params['bowtie2_align'], cutadapt.out.fastq, genome_ind )
}