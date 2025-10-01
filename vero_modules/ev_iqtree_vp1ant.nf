#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process iqtree_vp1{
    publishDir "${params.output}", mode: 'copy' // was commented - changed on 08/28/24
    input:
    //val samp    
    path vp1_aln
    //tuple  val(seq), path(R1), path(R1)
    //tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    path "iqtree_out/*.bionj"
    path "iqtree_out/*.ckp.gz"
    path "iqtree_out/*.iqtree"
    path "iqtree_out/*.log"
    path "iqtree_out/*.mldist"
    path "iqtree_out/*.model.gz"
    path "iqtree_out/*.treefile"
    //val samp
    //path("${seq}.skesa.fa")
    //path ("*")
    //E180-R2_S22/pilon_cons/E180-R2_S22_pilon.fasta
    //phylogenetic tree bwt given sample and public strains of EV
    script:
    """     
    # building tree from alignment
    mkdir -p iqtree_out
    iqtree -s ${vp1_aln} -pre iqtree_out/vp1 -af fasta -nt ${task.cpus}
       
    """
}

//${params.output}/clustalo_out/${vp1_aln}.fasta
// Size per contig (moodule loadseqkit)
// seqkit fx2tab --length --name --header-line file.fasta
//${params.output}/clustalo_out/${vp1_aln}.fasta