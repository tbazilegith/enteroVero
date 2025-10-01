process quast {
   input:
   val samp
   output:
   //path 'xfile.txt', emit: aLook
   val "${params.output}/${samp}"
   //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
    """     
    quast.py -o ${params.output}/${samp}/${samp}_assembly/quast_results/ ${params.output}/${samp}/pilon_cons/${samp}_pilon.fasta --min-contig 100
  
    """
}