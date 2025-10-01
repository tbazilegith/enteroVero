#!/usr/bin/env nextflow

// Processing reads to get metrics

process readsproc {
    input:
    val mypath
    path pyoutputs

    output:
    //stdout
    val mypath
    path pyoutputs

    """
    #echo ${mypath}
    genome=\$(cat ${pyoutputs} | cut -d "," -f 6) # genome in python asbl table, 6nd column
    #echo \${genome}
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    #echo \${samplename}
    
    run_assembly_shuffleReads.pl -gz ${mypath}/\${samplename}_clean_1.fq.gz ${mypath}/\${samplename}_clean_2.fq.gz > ${mypath}/\${samplename}_clean_shuffled.fq.gz
    run_assembly_readMetrics.pl ${mypath}/\${samplename}_clean_shuffled.fq.gz -e \${genome} > ${mypath}/\${samplename}_readMetrics.txt

    """
    }