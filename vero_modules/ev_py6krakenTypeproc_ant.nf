#!/usr/bin/env nextflow

// Processing kraken output for parsing EV gentotype


process kraken_gtypeproc {
    input:
    val mypath
    path pyoutputs

    output:
    //stdout
    val mypath
    path pyoutputs

    """
    #!/usr/bin/env python3

    import subprocess
    import pandas as pd
    import itertools
    import fileinput

    items = "${mypath}".strip().split("/")
    fa = open("${pyoutputs}", "a") # file to write report results
    filepath6 = "${mypath}/kraken_out/"+items[-1]+"_gtype.report"
    try:
        with open(filepath6, 'r') as kreport:
            #lines = kreport.readlines()
            for line in kreport:
                # rank code is the fourth field
                fields = line.strip().split('\t')
                if len(fields) >= 4 and fields[3].strip() == 'S2':
                    percentage = float(fields[0].strip())
                    taxonomic_name = fields[5].strip()
                         
                #l_parse = l.lstrip().rstrip().split("\t")
                #percent = l_parse[0]
                #tax_level = l_parse[3]
                #tax = l_parse[5].lstrip()
                #if tax_level == 'S1':
                    fa.write(","+str(taxonomic_name)+","+str(percentage))
                else:
                    fa.write(","+str(fields[5].strip())+","+str(float(fields[0].strip())))
                    
    except FileNotFoundError:
        percentage = "Check sample"
        taxonomic_name = "EVtype unclassified"
        fa.write(","+str(taxonomic_name)+","+str(percentage))
    except Exception as e:
        percentage = "Error occurred {e}"
        taxonomic_name = "Error occured"
        fa.write(","+str(taxonomic_name)+","+str(percentage))
        #f = open("${pyoutputs}", "a")
        #with open("${pyoutputs}", "a") as fa:
        #fa.write(","+str(tax)+","+str(percent))
    fa.close()
        
    """
}