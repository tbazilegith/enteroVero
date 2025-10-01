#!/usr/bin/env nextflow

// Indexing genome files for alignment purpose 
nextflow.enable.dsl = 2

process bwa_INDEX_each {
    //publishDir "${params.output}/bwa_index", mode: 'copy' // new changes   
    input:
    val samp
        
    output:
    val "${samp}"
               
    script:
    """
    bwa index -a is ${params.output}/${samp}/dnovo_out/${samp}_contigs.fa
       
    """
}