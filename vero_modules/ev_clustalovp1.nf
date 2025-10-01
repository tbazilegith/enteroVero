#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process clustalo_vp1{
    publishDir "${params.output}", mode: 'copy'
    //publishDir "${params.output}/clustalo_out", pattern: "*.fasta", mode: 'copy'
    //publishDir "${params.output}/clustalo_out", pattern: "*.fasta.mat", mode: 'copy'
    //publishDir "${params.output}/clustalo_out", pattern: "*.tree", mode: 'copy'

    input:
    path vp1_seqs    
    
    //tuple  val(seq), path(R1), path(R1)
    //tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    path "clustalo_out/*_aln.fasta"
    path "clustalo_out/*cat.fasta"
    path "clustalo_out/*.mat"
    path "clustalo_out/*.tree"
    //path "${params.output}/clustalo_out/*"
    //path params.output
    //val samp
    //path("${seq}.skesa.fa")
    //path("${params.output}/clustalo_out/*")

    //E180-R2_S22/pilon_cons/E180-R2_S22_pilon.fasta   

    script:
    """     
    # Multiple alignment btw sample-assemblies and ref genomes
    mkdir -p clustalo_out
        
    # concatenating all fasta into a multi-fasta
    cat ${vp1_seqs}/*.fasta > clustalo_out/vp1_cat.fasta
    
    # Aligning samples assemblies and vp1 
    clustalo -i clustalo_out/vp1_cat.fasta -o clustalo_out/vp1_aln.fasta --outfmt=fasta --distmat-out=clustalo_out/vp1_aln.fasta.mat --guidetree-out=clustalo_out/ev_vp1.tree --full
          
    """
}

// Size per contig (moodule loadseqkit)
// seqkit fx2tab --length --name --header-line file.fasta 