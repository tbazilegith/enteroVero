#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process mega_treeMod{
    //publishDir "${params.output}/", mode: 'copy'
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
    mkdir -p ${params.output}/${xp}/mega_out/mega_tree #
    #megacc -a ${params.mao2} -d ${params.output}/${xp}/mega_out/mega_aln/EV_vp1_aln.meg -o ${params.output}/${xp}/mega_out/mega_tree
    apptainer  exec docker://pegi3s/megax_cc:latest megacc -a ${params.mao4tree} -d ${params.output}/${xp}/mega_out/mega_aln/${xp}_aln.meg -o ${params.output}/${xp}/mega_out/mega_tree
    # May delete vp1 fasta sequences
    # rm ${params.output}/${xp}/ev_vp1/*.fasta 
       
    """
}

// Size per contig (moodule loadseqkit)
// seqkit fx2tab --length --name --header-line file.fasta
