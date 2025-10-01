#!/usr/bin/env nextflow

// Processing report for reads, speciesID, and Assembly

process reportproc {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    $/
    #!/usr/bin/env python3

    items = "${mypath}".strip().split("/")
    #parentpath = "/".join(items[:-1])
    
    with open("${pyoutputs}", "r") as aline:
        for line in aline:
            cells = line.rstrip().split(",")
	    #columns in pyoutputs table
            #results=[items[-1],cells[12],cells[11],cells[9],cells[7],cells[8],cells[10],cells[1],cells[2],cells[3],cells[4],cells[5],cells[6]] original table
            results = [items[-1],cells[12],cells[11],cells[9],cells[7],cells[8],cells[10],cells[1],cells[2],cells[3],cells[4],cells[5],cells[6],cells[13],cells[14]]
            #print(results)

    report = open("${mypath}"+"/report.txt", 'w')

    #Column headers
    #header = ['SampleID', 'SpeciesID_kraken', 'Kraken_percent', 'Num_clean_reads', 'Avg_readlength', 'Avg_read_qual', 'Est_coverage', 'Num_contigs', 'Longest_contig', 'N50', 'L50', 'Lotal	#_length', 'GC_content'] original table
    
    header = ['SampleID', 'SpeciesID_kraken', 'Kraken_percent', 'Num_clean_reads', 'Avg_readlength', 'Avg_read_qual', 'Est_coverage', 'Num_contigs', 'Longest_contig', 'N50', 'L50', 'Lotal_length', 'GC_content', 'KrEV_Type', 'Kr_%fragNtaxon']
    report.write("\t".join(header))
    report.write('\n')
    report.write('\t'.join(results))
    report.close()
    /$
}