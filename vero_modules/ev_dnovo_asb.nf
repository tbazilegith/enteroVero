#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.asbler = "megahit"

process dnovo_asbl {
    cpus = 12
    memory = 40.GB 
    //publishDir "${params.output}/skesa_dnovo", mode: 'copy'
    input:
    val samp    
    
    //tuple  val(seq), path(R1), path(R1)
    //tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    val samp
    //path("${seq}.skesa.fa")
    //path ("*")
    //errorStrategy 'ignore'
    //errorStrategy 'retry' // or 'ignore' to skip the failed task and continue
    //maxRetries 3 // Number of times to retry the task before failing

   errorStrategy {
        if(params.asbler == 'megahit') {
            'terminate' // megahit asbler, stop on error
        } else if(params.asbler == 'skesa') {
            'ignore'    // skesa, ignore errors
        } else{
	    'retry'
	}    
    }

    script:
    if( params.asbler == 'megahit' )
        """
        mkdir -p ${params.output}/${samp}/dnovo_out
        apptainer exec docker://biocontainers/megahit:1.2.9_cv1 megahit -1 ${params.output}/${samp}/${samp}_clean_1.fq.gz -2 ${params.output}/${samp}/${samp}_clean_2.fq.gz --out-prefix ${samp} -o ${params.output}/${samp}/megahit_out
        mv ${params.output}/${samp}/megahit_out/${samp}.contigs.fa ${params.output}/${samp}/dnovo_out/${samp}_contigs.fa
        """
    
    else if( params.asbler == 'skesa' )
        
        """     
        # Assembling reads
        mkdir -p ${params.output}/${samp}/dnovo_out
        apptainer exec docker://staphb/skesa skesa --reads ${params.output}/${samp}/${samp}_clean_1.fq.gz,${params.output}/${samp}/${samp}_clean_2.fq.gz > ${params.output}/${samp}/${samp}.skesa.fa
        mv ${params.output}/${samp}/${samp}.skesa.fa ${params.output}/${samp}/dnovo_out/${samp}_contigs.fa
        """

    else
        error "Invalid assembly tool: ${params.asbler}"
}

// Size per contig (moodule loadseqkit)
// seqkit fx2tab --length --name --header-line file.fasta
// core and memory ${task.cpus} ${task.memory}