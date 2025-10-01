#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process blast_nucl {
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
    # Blast of assembly
    mkdir -p ${params.output}/${samp}/blastn_results
    wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz -O - | tar -xz
    blastn -db ${params.blastnevdb}/${params.ev_ndbname} -query ${params.output}/${samp}/pilon_cons/${samp}_pilon.fasta -max_target_seqs 5 -outfmt "6 qseqid sseqid pident sscinames staxids evalue " -out ${params.output}/${samp}/blastn_results/${samp}_blastn.tsv
           
    """
}


