#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process bbtools {
    input:
    val samp

    output:
    //path 'xfile.txt', emit: aLook
    val "${samp}"
    //path "${params.output}/${samp}_trim_2.fastq", emit: trimR2
    
    script:

    """
    bbduk.sh in1=${params.output}/${samp}/${samp}_trim_1.fastq.gz in2=${params.output}/${samp}/${samp}_trim_2.fastq.gz out1=${params.output}/${samp}/${samp}_1.rmadpt.fq.gz out2=${params.output}/${samp}/${samp}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/${samp}/${samp}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo

    bbduk.sh in1=${params.output}/${samp}/${samp}_1.rmadpt.fq.gz in2=${params.output}/${samp}/${samp}_2.rmadpt.fq.gz out1=${params.output}/${samp}/${samp}_1.fq.gz out2=${params.output}/${samp}/${samp}_2.fq.gz outm=${params.output}/${samp}/${samp}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/${samp}/${samp}_phixstats.txt
   
    rm ${params.output}/${samp}/${samp}_trim*.fastq.gz
    rm ${params.output}/${samp}/${samp}*rmadpt.fq.gz
    """
}


// Adapter trimming, ktrim=r(right trimming or 3' adapters), k = kmer size, mink = look for shorter kmer length at the end of reads, hdist = hamming distance or allowed mismatch, ref = reference file, tbo trim based on overlap
// tpe trim both reads to the same length
//bbtools includes Illumina Truseq and Nextera adapters sequences in /bbmap/resources/adapters.fa
//kmer filtering - removing reads 31-mer match to PhiX, outm to catch reads matching a reference kmer, stats reports contaminant sequences