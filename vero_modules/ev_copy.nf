#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// Input is from pilon polished assembly and vp1 genes channel

process copy_files{
    input:
    val xp //samp
    output:
    val xp //samp
    
    script:
    """
    mkdir -p ${params.output}/${xp}/ev_vp1
    cp ${params.ev_ref_genomes}/*.fasta ${params.output}/${xp}/ev_vp1/
    cp ${params.output}/${xp}/pilon_cons/${xp}_pilon.fasta ${params.output}/${xp}/ev_vp1/
    """
}
