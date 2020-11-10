#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Module inclusions
include { fastqc } from '../../luslab-nf-modules/tools/fastqc/main.nf'
include { cutadapt } from '../../luslab-nf-modules/tools/cutadapt/main.nf'
include { bowtie2_align } from '../../luslab-nf-modules/tools/bowtie2/main.nf'

workflow qc_align {
    take: data
    take: genome_index
    take: module_params
    main:
        // Run fastqc
        fastqc( module_params.modules['fastqc'], data )

        // Adapter trimming
        cutadapt( module_params.modules['cutadapt'], data )

        // Align to genome
        bowtie2_align( module_params.modules['bowtie2_align'], cutadapt.out.fastq, genome_index )

    emit: bam = bowtie2_align.out.bam
    emit: fastqc_report = fastqc.out.report
    emit: cutadapt_report = cutadapt.out.report
    emit: bt2_report = bowtie2_align.out.report
}