#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process unmapped_proc {

    input:
        val samp

    output:
        val "${samp}",  emit: unmapped_outpath1

    script:
        """
        mkdir -p ${params.output}/${samp}/unmappings # to check
        samtools view -b -S ${params.output}/${samp}/mappings/${samp}.sam > ${params.output}/${samp}/unmappings/${samp}.bam
        samtools view --threads $task.cpus -b -F 2 ${params.output}/${samp}/unmappings/${samp}.bam > ${params.output}/${samp}/unmappings/${samp}_unmapped.bam
        samtools sort -n ${params.output}/${samp}/unmappings/${samp}_unmapped.bam -o ${params.output}/${samp}/unmappings/${samp}_umapped_sorted.bam
        samtools bam2fq -1 ${params.output}/${samp}/unmappings/${samp}_1.fq -2 ${params.output}/${samp}/unmappings/${samp}_2.fq -0 /dev/null -s /dev/null -n ${params.output}/${samp}/unmappings/${samp}_umapped_sorted.bam
        """
}

// /${samp}/mappings/${samp}.sam
//retrieve unmapped reads
//mkdir -p unmapped_bam
//unmappedR1.fastq and unmappedR2.fastq (both paired and unpaired)
//samtools view -b -f 4 bowtie_out/E175-R1_S10_L001_all.bam > unmapped_bam/E175-R1_S10_L001_unmapped.bam # good for single_end to

// or unmappedpairedR1.fastq (containing only paired R1 unmapped reads) and unmappedpairedR2.fastq (containing only paired R2 unmapped reads)
//samtools view --threads $PROCESSORS -b -F 2 in.bam > unmapped.bam (paired-end reads)
// sorting umapped reads
// samtools sort -n unmapped.bam -o umapped_sorted.bam

// bam2fasq  for unmapped reads (use sorted bam)
// mkdir -p unmapped_fastqs
// samtools bam2fq --no-aligned --force --strict -o unmapped_fastqs/E175-R1_S10_L001_unmapped#.fq bowtie_out/E175-R1_S10_L001_allsorted.bam
// samtools bam2fq -1 unmapped_fastqs/E175-R1_S10_L001_1.fq -2 unmapped_fastqs/E175-R1_S10_L001_2.fq -0 /dev/null -s /dev/null -n unmapped_bam/E175-R1_S10_L001_unmapped.bam
