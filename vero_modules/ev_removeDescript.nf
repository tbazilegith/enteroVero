#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// module to remove the description from the fasta file - and to keep only sequence ID
process rem_descr{
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
    mkdir -p ${params.output}/${xp}/mega_out #
        
    # concatenating all fasta (vp1 and sample assembly) into a multi-fasta
    cat ${params.output}/${xp}/ev_vp1/*.fasta > ${params.output}/${xp}/ev_vp1/${xp}_vp1_cat.fasta  
    remove_desc_arg.py --infasta ${params.output}/${xp}/ev_vp1/${xp}_vp1_cat.fasta  --outfasta ${params.output}/${xp}/ev_vp1/${xp}_vp1_corr.fasta    
   
    # May delete vp1 fasta sequences
    # rm ${params.output}/${xp}/ev_vp1/*.fasta 
       
    """
}

// Size per contig (moodule loadseqkit)
// seqkit fx2tab --length --name --header-line file.fasta
