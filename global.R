#Functions.R
#FUNCTIONS FOR TEXT INPUT
#require(tidyverse)f
require(shiny)
require(stringr)
require(doMC)
require(ggplot2)
require(plotly)
registerDoMC(85) #change to edit how many files %dopar% can run in parallel
require(shinyjs)
require(stringr)
require(tidyr)


#Read in frequency file for HLA alleles
hla <- read.delim('data/HLAtable.txt',header=T, stringsAsFactors=F)

#select from data/maps/GENE directory
select_maps <- list.files(path = paste("data/maps/"))

############################################################################
#TEXT INPUT PROTEIN
genetext_table <- data.frame(matrix(ncol = 3, nrow=1))
colnames(genetext_table) <- c("gene","unitprot","seq")

#read neo table
neo_seq_df <- read.csv("data/neomuts_table.csv")


#read total table of genes
#total_data <- read.csv("totaldata.csv", check.names = TRUE)

#validating text entries
# check_length <- function(input) {
#   aa_length <- str_length(input)
#   if (aa_length > 9) {
#     "aa length too long"
#   } else if (length < 9) {
#     "aa length too short"
#   } else {
#     NULL
#   }
#}






#Set Sequences
gene_seq_df <- read.csv("data/gene_seq.csv")
#Use MYCN Gene Frequency
genes <- c("TP53", "MYCN", "MLN", "BRAF", "PI3K")
gentab <- data.frame(matrix(ncol = 3, nrow=length(genes)))
colnames(gentab) <- c("gene","unitprot","seq")
gentab$gene <- genes


#CUSTOM NEOANTIGEN TABLE
peptab <- data.frame(matrix(ncol = 3, nrow=2))
colnames(peptab) <- c("gene","unitprot","seq")
peptab$gene <- c("WT1", "Mutant1") 

#LIBRARY NEOANTIGEN TABLE
neolibtab <- data.frame(matrix(ncol = 3, nrow=18))
colnames(neolibtab) <- c("gene","unitprot","seq")

#create peptide searchfile
createpep_searchfile <- function(sel_gene_df){
  sink(paste("data/NeoAntigens/wt_v_mut_netmhc.txt"))
  for (j in 1:nrow(sel_gene_df)){
    cat(paste(">", sel_gene_df$gene[j], "_", j, sep =""))
    cat("\n")
    cat(sel_gene_df$seq[j])
    cat("\n")
  }
  sink()
}

#insert text input pep seq into gene_seq_df
neotab <- data.frame(matrix(ncol = 3, nrow=nrow(neo_seq_df)))
colnames(neotab) <- c("gene","unitprot","seq")
neotab$gene <- neo_seq_df$gene
neotab$seq <- neo_seq_df$pep

#create LIBRARY NEOANTIGEN searchfile function
createneo_searchfile <- function(sel_neo_df){
  sink(paste("data/NeoAntigens/wt_v_mut_netmhc.txt"))
  for (i in 1:nrow(sel_neo_df)){
    for (j in 9:26) {
      cat(paste(">", sel_neo_df$gene[i], "_", sel_neo_df$substitution[i],
                "_", colnames(sel_neo_df)[j], sep =""))
      cat("\n")
      cat(paste(sel_neo_df[i, j]))
      cat("\n")
    }
  }
  sink()
} 


########################################################################
#protein/gene search
########################################################################
#creating searchfile for text input
createsearchfile <- function(sel_gene_df){
  for (j in 1:nrow(sel_gene_df)){
    a <- data.frame(matrix(ncol = 1))
    for (i in 1:(nchar(as.vector(sel_gene_df$seq[j]))-8)){
      a<-rbind(a,substr(sel_gene_df$seq[j], i, i+8))  
    }
    sink(paste("data/", sel_gene_df$gene[j], "netmhc.txt", sep=""))
    for (i in 2:nrow(a)){
      cat(paste(">",sel_gene_df$gene[j], "_", i-1,sep=""))
      cat("\n")
      cat(a$matrix.ncol...1.[i])
      cat("\n")
    }
    sink()
  }
}

#creating rules for
#A, C, D E F G H I K L M N P Q R S T V W Y and X (unknown)
allowed_aa <- c("A", "C", "D", "E", "F", "G", "H", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X")

#function to create binders
createbinders <- function(sel_gene_df, bindersfiles) {
  for (f in bindersfiles){
    a <- read.delim(f ,header=T,stringsAsFactors=F, sep="")
    b <- subset(a, a$BindLevel=="<=SB")
    f <-  f %>% str_replace(".*/", "")
    write.table(b, file=paste("data/", sel_gene_df$gene[1], "/peptides/binders/", f, sep=""), quote = F, row.names = F, sep ="\t")
  }
  print("Binder Creation Complete")
}


#function to create calculated HLA combined_table
createtab <- function(combined_table, sel_gene_df) {
  for (i in 1:ncol(combined_table)){
    r <- match(colnames(combined_table[i]), hla$Allele)
    combined_table["HLA_frequency", i] <- hla$HLA.Population.Frequency[r]
  }
  combined_table[is.na(combined_table)] <- as.numeric(0)
  
  for (i in 1:(nrow(combined_table)-1)){
    neofreqTCGA <- c()
    hlabinders <- c()
    
    for (j in 1:84){
      if (as.numeric(combined_table[i,j]) > 0){
        neofreqTCGA <- c(neofreqTCGA, combined_table["HLA_frequency", j])
        hlabinders <- c(hlabinders, substring(colnames(combined_table[j]), 5,9))
      }
      probTCGA <- 1
      if (length(neofreqTCGA) > 0){
        for (k in 1:length(neofreqTCGA)){
          probTCGA <- probTCGA*(1-as.numeric(neofreqTCGA[k]))
        }
      }
    }
    probTCGA <- 1-probTCGA
    combined_table$"HLA_frequency"[i] <- probTCGA
    combined_table$"Alleles_bound"[i] <- length(hlabinders)
    combined_table$"HLA_binders"[i] <- paste(unlist(hlabinders), collapse = ", ")
  }
  
  
  
  #code takes "_" and changes MYCN_1 to 1, need "xxxx_xxxx_wt1" to become "1"
  #code for creating table for normal tables (neo cust, gene cust/lib)
  # combined_table$position <- row.names(combined_table)
  # combined_table$gene <- gsub( "_.*$", "", row.names(combined_table))
  # for (i in 1:nrow(combined_table)){
  #   #column that produces amino acid position
  #   combined_table$position[i] <- sub('.*\\_', '', row.names(combined_table)[i])
  # }



  # for (i in 1:nrow(combined_table)){
  #   r <- match(gsub(" ", "", combined_table$gene[i]), sel_gene_df$gene)
  #   combined_table$aa[i] <- substr(sel_gene_df$seq[r], combined_table$position[i], combined_table$position[i])
  #   combined_table$pep[i] <- substr(sel_gene_df$seq[r], combined_table$position[i], (as.numeric(combined_table$position[i])+8))
  # }
  # write.csv(combined_table, paste("peptidemap.csv"))
  combined_table
  
}


