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
    items = "${mypath}".strip().split("/")
    filepath6="${mypath}/kraken_out/"+items[-1]+"_gtype.report"
    with open(filepath6, 'r') as kreport:
        lines = kreport.readlines()
        for l in lines:
            l_parse = l.lstrip().rstrip().split("\t")
            percent = l_parse[0]
            tax_level = l_parse[3]
            tax = l_parse[5].lstrip()
            if tax_level == 'S1':
                break
        #f = open("${pyoutputs}", "a")
    with open("${pyoutputs}", "a") as fa:
        fa.write(","+str(tax)+","+str(percent))
        #f.close()
        
    """
}