#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// Taxonomic classification using k-mers

process kraken {
    input:
        val samp
    output:
        //stdout
        val '${samp}'
    script:    
    """   
    mkdir -p ${params.output}/${samp}/kraken_out/
    #kraken2 --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq

    kraken2 --db ${params.db} --use-names --threads $task.cpus --report ${params.output}/${samp}/kraken_out/${samp}.report --output ${params.output}/${samp}/kraken_out/${samp}_kraken.out --paired ${params.output}/${samp}/${samp}_clean_1.fq.gz ${params.output}/${samp}/${samp}_clean_2.fq.gz
    """
}