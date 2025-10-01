#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Converting dna to protein

process dna2protein_proc {
    input:
    val mypath
    //path pyoutputs
    output:
    //stdout
    val mypath
    //path pyoutputs
    script:
    """
    #items = "${mypath}".strip().split("/")
    #filepath3="${mypath}/pilon_cons/${mypath}_pilon.fasta"
    dna2proteins_v2.py --fastapath ${params.output}/${mypath}/pilon_cons/${mypath}_pilon.fasta --outpro ${params.output}/${mypath}/pilon_cons/${mypath}_protein.faa 
    
    """
}

