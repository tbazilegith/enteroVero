#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process scrub {
    //publishDir "${params.output}/scrub_out", mode: 'copy'

    input:
        val seq    
    
    output:
    //path("*")
        //val "${params.output}/${seq}"
        val "${seq}"
    script:
    """
    # decompress
    gzip -d ${params.output}/${seq}/${seq}_1.fastq.gz 
    gzip -d ${params.output}/${seq}/${seq}_2.fastq.gz      
    # interleaving - combine reads
    #reformat.sh in1=${params.output}/${seq}/${seq}_1.fq in2=${params.output}/${seq}/${seq}_2.fq out=${params.output}/${seq}/${seq}_interleaved.fq

    # scrubbing pe reads
    scrub.sh -i ${params.output}/${seq}/${seq}_1.fastq -o ${params.output}/${seq}/${seq}_1_clean.fastq
    scrub.sh -i ${params.output}/${seq}/${seq}_2.fastq -o ${params.output}/${seq}/${seq}_2_clean.fastq
    # deinterleave
    #reformat.sh in=${params.output}/${seq}/${seq}.fq.clean out1={params.output}/${seq}/${seq}_1.fq.clean out2=${params.output}/${seq}/${seq}_2.fq.clean
    # compress
    gzip ${params.output}/${seq}/${seq}_1_clean.fastq
    gzip ${params.output}/${seq}/${seq}_2_clean.fastq
    """
}