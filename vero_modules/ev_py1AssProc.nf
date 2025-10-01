#!/usr/bin/env nextflow

// processing assembly statistics

process pyAsbProc {
    input:
    val mypath

    output:
    //stdout
    val "${mypath}", emit: outputpath2
    path "pyoutputs.txt", emit: pyoutputs
    $/
    #!/usr/bin/env python3
    import sys
    import os
    import subprocess
    import argparse
    import datetime
    import fnmatch
    import re
    import pandas as pd

    items = "${mypath}".strip().split("/")
    
    # Accessing the quast results
    filepath1 = "${mypath}"+"/"+items[-1]+"_assembly/quast_results/report.tsv"
    # Open the quast results in the df
    df = pd.read_table(filepath1, sep="\t")
    assem = list(df.columns)[1]
    #print(assem)
    contigs = df[assem][12].astype(int)
    long_contig = df[assem][13].astype(int)
    n50 = df[assem][16].astype(int)
    l50 = df[assem][18].astype(int)
    genome = df[assem][14].astype(int)
    gc = df[assem][15].astype(int)
    #print(genome)

    # File for writing output
    f = open("pyoutputs.txt", "w")
    f.write(str(assem)+","+str(contigs)+","+str(long_contig)+","+str(n50)+","+str(l50)+","+str(genome)+","+str(gc))
    f.close

    /$
    }