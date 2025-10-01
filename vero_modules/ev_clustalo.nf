#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process clustalo{
    //publishDir "${params.output}/clustalo_out", mode: 'copy'
    input:
    val xp //samp    
    
    //tuple  val(seq), path(R1), path(R1)
    //tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    val xp //samp
    
    //path("${seq}.skesa.fa")
    //path ("*")
    //E180-R2_S22/pilon_cons/E180-R2_S22_pilon.fasta   

    script:
    """     
    # Multiple alignment btw sample-assemblies and ref genomes
    mkdir -p ${params.output}/${xp}/clustalo_out #
        
    # concatenating all fasta into a multi-fasta
    #cat ${params.output}/${xp}/ev_vp1/*.fasta > ${params.output}/${xp}/ev_vp1/${xp}_vp1_cat.fasta  
    
    clustalo -i ${params.output}/${xp}/ev_vp1/${xp}_vp1_cat.fasta -o ${params.output}/${xp}/clustalo_out/${xp}_evpub.fasta --outfmt=fasta --distmat-out=${params.output}/${xp}/clustalo_out/${xp}_evpub.fasta.mat --guidetree-out=${params.output}/${xp}/clustalo_out/${xp}_vp1.tree --full --force --threads ${task.cpus}
    # May delete vp1 fasta sequences
    # rm ${params.output}/${xp}/ev_vp1/*.fasta 
       
    """
}

// Size per contig (moodule loadseqkit)
// seqkit fx2tab --length --name --header-line file.fasta
