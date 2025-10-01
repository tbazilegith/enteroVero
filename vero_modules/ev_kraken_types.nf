#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Kraken Taxonomic classification using k-mers
// Classification of sample by EV genotype

process kraken4gtype {
    input:
        val samp
    output:
        //stdout
        val '${samp}'
    script:    
    """   
    #mkdir -p ${params.output}/${samp}/kraken_out/
    #kraken2 --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq

    kraken2 --db ${params.vp1db} --use-names --threads $task.cpus --report ${params.output}/${samp}/kraken_out/${samp}_gtype.report --output ${params.output}/${samp}/kraken_out/${samp}_gtype_kraken.out ${params.output}/${samp}/pilon_cons/${samp}_pilon.fasta
    """
}