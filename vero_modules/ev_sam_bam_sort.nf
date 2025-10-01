#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process samtools{
    //publishDir "${params.output}/samtools_sortedbam", mode: 'copy'   
    input:
    val samp 
    
    output:
    val samp       

    script:
    """
    samtools view -b ${params.output}/${samp}/mappings/${samp}.sam | samtools sort - > ${params.output}/${samp}/mappings/${samp}.sorted.bam
    samtools index ${params.output}/${samp}/mappings/${samp}.sorted.bam
    """
    }