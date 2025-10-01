#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
==================================================================================
ev_gen
A Workflow to identify genotypes of Enteroviruses
==================================================================================
https://github.com/
#### Author
Tassy J. Bazile <tassy.bazile@flhealth.gov> - https://github.com/

--------------------------------------------------------------------------------
*/

/*
==================================================================================
Help Message
==================================================================================
*/

def helpMessage() {
    log.info"""
    >>
Usage:

The typical command for running the pipeline is as follows:
nextflow run ev_gen.nf --input_reads '*_{1,2}.fastq.gz' --outdir 'output path'

--input_reads   Path to short-read fastq files
--outdir        Chosen name of the output directory
<<
""".stripIndent()
}

//Show help message

if (params.help){
    helpMessage()
        exit 0
	}
/*
==================================================================================
Parameters
==================================================================================
*/
//params.input_reads = "" // format to provide paired-end reads  *_{1,2}.fastq') //params.genome = ""
//params.outdir = "$baseDir/fastqc_results"




/*
==================================================================================
Channels
==================================================================================
*/

//reads_ch = Channel.fromFilePairs(params.input_reads, checkIfExists: true )


/*
==================================================================================
Processes
==================================================================================
*/

process FASTQC {

    input:
    val samp

    output:
    val "${samp}",  emit: fastqc_outpath1

    script:
    """
    mkdir -p ${params.output}/${samp}
    cp ${params.input}/${samp}_*.fastq.gz ${params.output}/${samp}
    fastqc ${params.output}/${samp}/${samp}_1.fastq.gz ${params.output}/${samp}/${samp}_2.fastq.gz
    # Adding fastqc output in the same dir as input sample 
    mv ${params.output}/${samp}/${samp}_1_fastqc.html ${params.output}/${samp}/${samp}_1_original_fastqc.html
    mv ${params.output}/${samp}/${samp}_1_fastqc.zip ${params.output}/${samp}/${samp}_1_original_fastqc.zip
    mv ${params.output}/${samp}/${samp}_2_fastqc.html ${params.output}/${samp}/${samp}_2_original_fastqc.html
    mv ${params.output}/${samp}/${samp}_2_fastqc.zip ${params.output}/${samp}/${samp}_2_original_fastqc.zip
    """
}



/*
==================================================================================

Workflow
==================================================================================


workflow {

    FASTQC( reads_ch )
}
*/