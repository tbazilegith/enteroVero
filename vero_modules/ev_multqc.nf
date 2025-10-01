//E128_S9_clean_2_fastqc.zip
// Aggregating qc results accross samples into a single report
process multiqc {
    input:
    val samp
    output:
    //path 'xfile.txt', emit: aLook
    val "${samp}"
    //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
    script:
    """     
   
    multiqc ${params.output}/${samp}/${samp}_clean_*_fastqc.zip -o ${params.output}/${samp}
     
    """
    }