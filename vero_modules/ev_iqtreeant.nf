#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process iqtree{
    //publishDir "${params.output}/${samp}/iqtree_out/", mode: 'copy' // was commented - changed on 08/28/24
    input:
    val samp    
    
    //tuple  val(seq), path(R1), path(R1)
    //tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    val samp
    //path("${seq}.skesa.fa")
    //path ("*")
    //E180-R2_S22/pilon_cons/E180-R2_S22_pilon.fasta
    //phylogenetic tree bwt given sample and public strains of EV
    script:
    """     
    # building tree from alignment
    mkdir -p ${params.output}/${samp}/iqtree_out/
    iqtree -s ${params.output}/${samp}/clustalo_out/${samp}_evpub.fasta -pre ${params.output}/${samp}/iqtree_out/${samp}_evpub -af fasta -nt ${task.cpus}
        
    """
}

// Size per contig (moodule loadseqkit)
// seqkit fx2tab --length --name --header-line file.fasta 