#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
//params.mypath = ""
//params.vp1_metdt = ""
//params.output = "my_temp_results"


process treeAnnotTemp {
    //publishDir "${params.output}", mode: 'copy'
    input:
    val mypath

    output:
    val mypath

    script:
    """
        
    #!/usr/bin/env Rscript
    
    library("ape")
    library(ggtree)
    library("ggplot2")
    #library(writexl)
    library(readxl)
    library(RColorBrewer)
    library(ggnewscale)
    library(tidyverse)
    library(fs)# to find files in path
    library(TAF) # to make dir
    #items = "${mypath}".strip().split("/")
    
    # ouput directory
    #mkdir("tree_annot")
    
    # Accessing the file path
    #split_path  <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))
    #lsubdir <- split_path("${mypath}")
    #samp_name <- lsubdir[1]

    phyl_tree <- list.files(path = "${params.output}/${mypath}", pattern = ".treefile|.nwk", include.dirs= TRUE, full.names=TRUE, recursive = TRUE)

    # function iq_meg
    iq_meg <- function(xtree,dir4metdt){
        tree01 <- read.tree(xtree)
     
        # Getting tree tip labels *****
        vp1_accession <- tree01[4] # TipLabels(tree01)
        # Ordering
        vp1_accession2 <- vp1_accession[[1]] # as vector of characters , sort(as.vector(vp1_accession))
    
        # Arranging tree tips into df
        #vp1_accession3 <- as.vector(vp1_accession2)
        df_treetip <- as.data.frame(vp1_accession2)
    
        # Extracting, filter sample Name/contigs from vp1_accession vector or treetip labels
        df_ext <- df_treetip %>% filter(str_detect(vp1_accession2, 'Contig_'))
  
        #safe df as a csv file
        #write.csv(df_ext,"${params.output}/${mypath}/EV_vp1metadata.csv", row.names = FALSE)
        
        # Adding metadata columns to df_ext
        namevector <- c("AA_Accession","Species_N_Type", "Genotype", "Species")
        df_ext[,namevector] <- NA
    
        # Changing the name of colum  "vp1_accession3" to Nucleotide_Accession in df_ext
        df_ext <- df_ext %>% relocate(Nucleotide_Accession = vp1_accession2)
    
        # Read VP1 metadata df
        #df_evtypes <- read_excel("${params.output}/EV_RV_vp1GenProt.xlsx", na = "NA")
        df_evtypes <- read_excel("${params.vp1_metdt}", na = "NA")
        df_evtypes <- as.data.frame(df_evtypes) # write as a df
        df_evtypes[df_evtypes == ''] <- "NA" # dealing with NA
    
        # First column of df should  match the external node/tip labels of the tree.
        df_evtypes2 <- df_evtypes%>%relocate(Nucleotide_Accession)

        # Relocate colums AA_Accession next to Nucleotide_Accession
        df_evtypes2 <- df_evtypes2 %>% relocate(AA_Accession, .after = Nucleotide_Accession)
    
        #Renaming colunm Species-type
    
        #colnames(df_evtypes2)[colnames(df_evtypes2) == 'Species-type'] <- "Species_Type"
    
        #Duplicate the nucleotide accession column(has all sequence IDs) to show sampleID in clade species types

        df_ext1 = df_ext # Original treetip labels with only contigs of sample
    
        df_ext1 <- df_ext1 %>% mutate(Sample_Nuacc = Nucleotide_Accession)
    
        # Select needed columns (Sample_Nuacc replaces Species_type)
        df_ext2 <- select(df_ext1, Nucleotide_Accession, AA_Accession, Sample_Nuacc, Genotype, Species)
    
        #Change column Species_Nuacc name into Species_Type
        colnames(df_ext2)[colnames(df_ext2) == "Sample_Nuacc"] <- "Species_N_Type"
    
        # Combining df_ext2 and ddf_evtypes2
        df_treeNdata2 <- rbind(df_ext2, df_evtypes2)
    
        # Changing tree tips based on df Species_Type column to see clade of sample types
        tree02 <- tree01
        #treetl <- TipLabels(tree02)
        tree02[4][[1]] <- df_treeNdata2[[3]][match(tree02[4][[1]], df_treeNdata2[[1]])] # colunm 3(new tree tips) and colum 1 (original tree tips)

        # Annotated tree
        atree <- ggtree(tree02) %<+% df_treeNdata2 + geom_tiplab(align=TRUE, linesize=.3, size=1.7) + theme_tree2() + xlim(NA, 12) + ggtitle("Sample ${mypath} and VP1 trains of EV and RV")
        
        # Save metadata files
        file_name_1 <- paste0(dir4metdt, "/${mypath}_contigNvp1.csv")
        file_name_2<- paste0(dir4metdt, "/EV_vp1metadata.csv")
        write.csv(df_treeNdata2, file = file_name_1, row.names = FALSE)
        write.csv(df_evtypes2, file = file_name_2, row.names = FALSE)

        return(atree)
    }

    dir4metdt <- "${params.output}/${mypath}"
    for (i in 1:length(phyl_tree)){
        if (endsWith(phyl_tree[[i]], ".treefile")){
        pp = iq_meg(phyl_tree[[i]],dir4metdt=dir4metdt)
        mkdir("${params.output}/${mypath}/iq_treeAnnot")
        ggsave("${params.output}/${mypath}/iq_treeAnnot/${mypath}_iqtree.pdf", width = 8.5, height = 11, units = "in")
        }else{
        qq = iq_meg(phyl_tree[[i]],dir4metdt=dir4metdt)
        mkdir("${params.output}/${mypath}/mega_treeAnnot")
        ggsave("${params.output}/${mypath}/mega_treeAnnot/${mypath}_megatree.pdf", width = 8.5, height = 11, units = "in")
        }
    }

    
    
    
    #write.csv(df_treeNdata2,"${params.output}/EV_vp1NsampleMdt.csv", row.names = FALSE)
    """
}

/*
workflow {
}
comments     
*/


