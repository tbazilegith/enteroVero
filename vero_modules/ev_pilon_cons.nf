#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process pilon {
    //publishDir "${params.output}/skesa_dnovo", mode: 'copy'
    input:
    val samp    
    
    //tuple  val(seq), path(R1), path(R1)
    //tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    //path("${seq}.skesa.fa")
    //path ("*")
    val samp
    script:
    """     
    # Assembling reads
    mkdir -p ${params.output}/${samp}/pilon_cons

    java -jar -Xmx16G /pilon/pilon.jar --genome ${params.output}/${samp}/dnovo_out/${samp}_contigs.fa --changes --frags ${params.output}/${samp}/mappings/${samp}.sorted.bam --output ${samp}_pilon --threads ${task.cpus} --outdir ${params.output}/${samp}/pilon_cons/     
    """
}