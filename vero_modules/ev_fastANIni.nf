#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process fast_ani1{
    //publishDir "${params.output}/skesa_dnovo", mode: 'copy'
    input:
    val samp    
    
    //tuple  val(seq), path(R1), path(R1)
    //tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    val samp
    //path("${seq}.skesa.fa")
    //path ("*")
    //E180-R2_S22/pilon_cons/E180-R2_S22_pilon.fasta   
    script:
    """     
    # Assembling reads
    mkdir -p ${params.output}/${samp}/fast_ani/
   
    ls ${params.ev_ref_genomes}/*.fasta > ${params.output}/${samp}/ev_strains.txt
    fastANI -q ${params.output}/${samp}/pilon_cons/${samp}_pilon.fasta --rl ${params.output}/${samp}/ev_strains.txt -o ${params.output}/${samp}/${samp}.ani.out --threads ${task.cpus}     
   
    mv ${params.output}/${samp}/${samp}.ani.out ${params.output}/${samp}/fast_ani/
    
    """
}

// Size per contig (moodule loadseqkit)
// seqkit fx2tab --length --name --header-line file.fasta 