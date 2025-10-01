#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process trimmomatic {
    input:
    val samp
    output:
    //path 'sampfile.txt', emit: aLook
    val "${samp}"
   //path "${params.output}/${samp}_trim_2.fastq", emit: trimR2
    script:
    """     
    trimmomatic PE -phred33 -trimlog ${params.output}/${samp}/${samp}.log ${params.output}/${samp}/${samp}_1_clean.fastq.gz ${params.output}/${samp}/${samp}_2_clean.fastq.gz ${params.output}/${samp}/${samp}_trim_1.fastq.gz ${params.output}/${samp}/${samp}_unpaired_trim_1.fastq.gz ${params.output}/${samp}/${samp}_trim_2.fastq.gz ${params.output}/${samp}/${samp}_unpaired_trim_2.fastq.gz ILLUMINACLIP:/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:71 TRAILING:20 > ${params.output}/${samp}/${samp}_trimstats.txt   
    rm ${params.output}/${samp}/${samp}_unpaired_trim_*.fastq.gz
    rm ${params.output}/${samp}/${samp}_1.fastq ${params.output}/${samp}/${samp}_2.fastq
    """
}