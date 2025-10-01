#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process scrub2 {
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
    gzip -d ${params.output}/${seq}/${seq}_1.fq.gz 
    gzip -d ${params.output}/${seq}/${seq}_2.fq.gz      
    # interleaving - combine reads
    #reformat.sh in1=${params.output}/${seq}/${seq}_1.fq in2=${params.output}/${seq}/${seq}_2.fq out=${params.output}/${seq}/${seq}_interleaved.fq

    # scrubbing pe reads
    scrub.sh -i ${params.output}/${seq}/${seq}_1.fq -o ${params.output}/${seq}/${seq}_clean_1.fq
    scrub.sh -i ${params.output}/${seq}/${seq}_2.fq -o ${params.output}/${seq}/${seq}_clean_2.fq
    # deinterleave
    #reformat.sh in=${params.output}/${seq}/${seq}.fq.clean out1={params.output}/${seq}/${seq}_1.fq.clean out2=${params.output}/${seq}/${seq}_2.fq.clean
    # compress
    gzip ${params.output}/${seq}/${seq}_clean_1.fq
    gzip ${params.output}/${seq}/${seq}_clean_2.fq
    """
}