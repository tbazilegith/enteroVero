#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
==================================================================================
EnteroVero
A Workflow to identify genotypes of Enteroviruses
==================================================================================
https://github.com/
#### Author
Tassy J. Bazile <tassy.bazile@flhealth.gov> - https://github.com/BPHL-Molecular/enteroVero

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
nextflow run enterovero_wf.nf --input_reads '*_{1,2}.fastq.gz' --outdir 'output path'

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

//Step1:input data files

def L001R1Lst = []
def sampleNames = []
myDir = file("$params.input")

myDir.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}
L001R1Lst.sort()
L001R1Lst.each{
    def samp = it.minus("_1.fastq.gz")
        //println samp
    sampleNames.add(samp)
}
//println L001R1Lst
//println sampleNames


//Step2: process the inputed data
A = Channel.fromList(sampleNames)


/*
==================================================================================
Channels
==================================================================================
*/

//reads_ch = Channel.fromFilePairs(params.input_reads, checkIfExists: true )


/*
==================================================================================
Module Importation
==================================================================================
*/

include { FASTQC as FASTQC_ev } from './vero_modules/ev_gen_fastqc.nf'
include { trimmomatic as trimmomatic_ev } from './vero_modules/ev_trim.nf'
include { bbtools as bbtools_ev } from './vero_modules/ev_bbtools.nf' 
include {scrub as scrub_ev } from './vero_modules/ev_scrub.nf'
include {FASTQC2 as FASTQC2_clean } from './vero_modules/ev_gen_fastqc2.nf'
include {scrub2} from './vero_modules/ev_scrub2.nf'
include {kraken as kraken_ev} from './vero_modules/ev_kraken.nf'
include {dnovo_asbl as dnovo_asbl_ev} from './vero_modules/ev_dnovo_asb.nf'
include {bwa_INDEX_each} from './vero_modules/ev_bwaIndexEach.nf'
include {bwa_aln_proc} from './vero_modules/ev_aln_proc.nf'
include {samtools} from './vero_modules/ev_sam_bam_sort.nf'
include {pilon} from './vero_modules/ev_pilon_cons.nf'
include {quast} from './vero_modules/ev_quast.nf' // added on 04/23
include {pyAsbProc} from './vero_modules/ev_py1AssProc.nf' // added on 04/23
include {readsproc} from './vero_modules/ev_readsproc.nf' // added on 04/23
include {multiqc} from './vero_modules/ev_multqc.nf' //added on 04/23
include {readstat} from './vero_modules/ev_py2readStat.nf'
include {krakenproc} from './vero_modules/ev_py3Krakenproc.nf'
include {reportproc} from './vero_modules/ev_py4ReportProc.nf' //commented on 08/30/2024
include {unmapped_proc} from './vero_modules/ev_unmapped.nf'
include {fast_ani1} from './vero_modules/ev_fastANIni.nf'
include {copy_files} from './vero_modules/ev_copy.nf' // added on 01/15/2025
include {rem_descr} from './vero_modules/ev_removeDescript.nf' // moved on 02/25/25
include {clustalo} from './vero_modules/ev_clustalo.nf'
include {iqtree} from './vero_modules/ev_iqtree.nf'
include {dna2protein_proc} from './vero_modules/ev_dna2protein2.nf'
include {blast_nucl} from './vero_modules/ev_blast_n.nf'
include {blast_prot} from './vero_modules/ev_blast_p.nf'
include {treeAnnotProc} from './vero_modules/ev_phytree_annot2.nf'
include {clustalo_vp1 } from './vero_modules/ev_clustalovp1.nf'
include {iqtree_vp1 } from './vero_modules/ev_iqtree_vp1.nf'
include {treeVP1Annot} from './vero_modules/ev_phytree_VP1annot.nf'
include {mega_aln} from './vero_modules/ev_mega_aln.nf'
include {mega_treeMod} from './vero_modules/ev_mega_tree.nf'
include {kraken4gtype} from './vero_modules/ev_kraken_types.nf'
include {kraken_gtypeproc} from './vero_modules/ev_py6krakenTypeproc.nf'
// Reuse Process
include {treeAnnotTemp} from './vero_modules/ev_treeAnnot_mod12.nf'
//include {treeAnnotTemp as treeAnnotTemp01} from './ev_treeAnnot_mod1.nf'
//include {treeAnnotTemp as treeAnnotMega} from './ev_treeAnnot_mod.nf'

/*
==================================================================================

Workflow
==================================================================================
*/

workflow {
// Quality control
ch_fastqc = FASTQC_ev(A)
ch_scrub = scrub_ev(ch_fastqc)
ch_scrub.view()
ch_trim = trimmomatic_ev(ch_scrub)
ch_bbt = bbtools_ev(ch_trim)
ch_scrub2 = scrub2(ch_bbt)//added on 24/04/25
ch_fast_clean = FASTQC2_clean(ch_scrub2)
ch_multqc = multiqc(ch_fast_clean) // all qc

// Taxonomic Classification
ch_kraken =  kraken_ev(ch_fast_clean)

// Assembly and consensus
ch_dnovo = dnovo_asbl_ev(ch_fast_clean)
ch_bwaidx_each = bwa_INDEX_each(ch_dnovo)
ch_join = ch_fast_clean.join(ch_bwaidx_each) // sequences  and  genome index

//Alignment
ch_aln = bwa_aln_proc(ch_join) // bwa alignment
ch_sam_bam = samtools(ch_aln) // sam to  bam and index
ch_join_cons = ch_sam_bam.join(ch_dnovo)

ch_pilon = pilon(ch_join_cons)
ch_quast = quast(ch_pilon) // assembly metrics
ch_ablproc = pyAsbProc(ch_quast)// parsing asb metrics in a table file - table file created
ch_readproc = readsproc(ch_ablproc[0],ch_ablproc[1]) //path and table file(columns appended)
ch_readst =  readstat(ch_readproc[0], ch_readproc[1])// parsing (appending)reads metrics to table file ([0]path, [1] table file) 

ch_fastANIni = fast_ani1(ch_pilon) // added on 08/08/24
ch_kreport = krakenproc(ch_readst[0], ch_readst[1])// parsing(appending) Kraken report in table file

ch_unmap = unmapped_proc(ch_aln) //new

ch_copy_gen = copy_files(ch_pilon) // copy sample assemblies and VP1

// Remove long description line from sequence ID
ch_rem_descr = rem_descr(ch_copy_gen)

// Multiple sequence alignment
ch_clustalo = clustalo(ch_rem_descr)//02/25/2025
//ch_clustalo = clustalo(ch_copy_gen)

ch_iqtree = iqtree(ch_clustalo) // tree ev-type vp1 and samples
ch_blastn = blast_nucl(ch_pilon)
ch_dna2protein = dna2protein_proc(ch_pilon)
ch_blastp = blast_prot(ch_dna2protein)
ch_treeAnnot = treeAnnotProc(ch_iqtree)

// tree with all vp1
vp1_ch = Channel.fromPath(params.ev_ref_genomes)
ch_clust_vp1 = clustalo_vp1(vp1_ch)
ch_iqtree_vp1 = iqtree_vp1(ch_clust_vp1[0]) // tree ev-type vp1 only
ch_vp1_tree = treeVP1Annot(ch_iqtree_vp1[6]) // check 02/05

// Mega
ch_mega_aln = mega_aln(ch_rem_descr) //02/25/2025 Remove long description line from sequence ID 
ch_mega_tree = mega_treeMod(ch_mega_aln)

// Genotyping with Kraken
ch_kraken_type = kraken4gtype(ch_pilon) // classify sample by genotypes
ch_krtypeproc = kraken_gtypeproc(ch_kreport[0],ch_kreport[1]) // Append type column to previous table file
// Sum report
ch_finreport = reportproc(ch_krtypeproc[0], ch_krtypeproc[1]) // New final report 03/06/25

// Reuse Process

//ch_annot = treeAnnotTemp(ch_iqtree) // tree ev-type vp1 and samples
ch_annot2 = treeAnnotTemp(ch_mega_tree)

}


