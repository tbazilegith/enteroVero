#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process bwa_aln_proc {
    //publishDir "${params.output}/bwa_align", mode: 'copy'
    input:
    
    val samp    
    
    output:
    val samp
    
    script:
    """
    mkdir ${params.output}/${samp}/mappings     
    # Aligning pe reads
    bwa mem -t $task.cpus ${params.output}/${samp}/dnovo_out/${samp}_contigs.fa ${params.output}/${samp}/${samp}_clean_1.fq.gz ${params.output}/${samp}/${samp}_clean_2.fq.gz > ${params.output}/${samp}/mappings/${samp}.sam
        
    """
    }