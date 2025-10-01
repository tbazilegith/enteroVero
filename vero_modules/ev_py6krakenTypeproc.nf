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
    s2_data =[]
    u_data = []     
    filepath6 = "${mypath}/kraken_out/"+items[-1]+"_gtype.report"
    try:
        with open(filepath6, 'r') as kreport:
            #lines = kreport.readlines()
            for line in kreport:
                # rank code is the fourth field
                fields = line.strip().split('\t')
                #if len(fields) >= 4 and fields[3].strip() == 'S2':
                if len(fields) >= 6:
                    percentage = float(fields[0].strip())
                    taxonomic_level = fields[3].strip()
                    taxonomic_name = fields[5].strip()
                    if taxonomic_level == 'S2':
                        s2_data.append(f"{percentage}\t{taxonomic_name}")
                        #fa.write(","+str(taxonomic_name)+","+str(percentage))
                        break
                    elif taxonomic_level == 'U':
                        u_data.append(f"{percentage}\t{taxonomic_name}")
                        #fa.write(","+str(taxonomic_name)+","+str(percentage))
        if s2_data:
            fa.write(","+str(taxonomic_name)+","+str(percentage))
      
        elif u_data:
            fa.write(","+str(taxonomic_name)+","+str(percentage))   
        else:
            percentage = None
            taxonomic_name = None
            fa.write(","+str(taxonomic_name)+","+str(percentage)) 
              
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
        
    fa.close()
        
    """
}