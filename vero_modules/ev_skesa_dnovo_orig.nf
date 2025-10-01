#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process skesa_asbl {
    //publishDir "${params.output}/skesa_dnovo", mode: 'copy'
    input:
    val samp    
    
    //tuple  val(seq), path(R1), path(R1)
    //tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    val samp
    //path("${seq}.skesa.fa")
    //path ("*")
   
    script:
    """     
    # Assembling reads
    mkdir -p ${params.output}/${samp}/skesa_dnovo
    skesa --reads ${params.output}/${samp}/${samp}_clean_1.fq.gz,${params.output}/${samp}/${samp}_clean_2.fq.gz --cores ${task.cpus} --memory ${task.memory} > ${params.output}/${samp}/${samp}.skesa.fa
    mv ${params.output}/${samp}/${samp}.skesa.fa ${params.output}/${samp}/skesa_dnovo/
    
    """
}

// Size per contig (moodule loadseqkit)
// seqkit fx2tab --length --name --header-line file.fasta 