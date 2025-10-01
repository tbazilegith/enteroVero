#!/usr/bin/env nextflow

// Processing kraken output for parsing


process krakenproc {
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
    filepath3="${mypath}/kraken_out/"+items[-1]+".report"
    with open(filepath3, 'r') as kreport:
        lines = kreport.readlines()
        for l in lines:
            l_parse = l.lstrip().rstrip().split("\t")
            percent = l_parse[0]
            tax_level = l_parse[3]
            tax = l_parse[5].lstrip()
            if tax_level == 'S':
                break
    f = open("${pyoutputs}", "a")
    f.write(","+str(percent)+","+str(tax))
    f.close
    """
}