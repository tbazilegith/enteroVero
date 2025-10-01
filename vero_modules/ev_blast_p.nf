#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process blast_prot{
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
    # Assembling reads ; wget is to download and inflate taxdb.tar.gz file to having scientific name 
    mkdir -p ${params.output}/${samp}/blastp_results
    wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz -O - | tar -xz
    blastp -db ${params.blastpevdb}/${params.ev_pdbname} -query ${params.output}/${samp}/pilon_cons/${samp}_protein.faa -max_target_seqs 5 -outfmt "6 qseqid sseqid pident sscinames staxids evalue" -out ${params.output}/${samp}/blastp_results/${samp}_blastp.tsv          
    """
}