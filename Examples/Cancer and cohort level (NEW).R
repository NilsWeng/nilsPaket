#Cancer and cohort level 
rm(list=ls())
gc()
#Library
library(nilsPaket)
library(dplyr)
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)



#Set wd to folder containing all data
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)


#Read cosmic signatures
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))
#cosmic_signatures <- as.matrix(read.table("cosmic_signatures.txt",header=TRUE))


#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output")
#mutational_matrix <- read.table("COMPLETE_MUT_MAT.txt",header=TRUE)


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output/NEW_COHORT")
mut_mat_files <- list.files()
mut_mat_files <- mut_mat_files[grep("_mut_matrix.txt",mut_mat_files)]




#mut_mat_files <- mut_mat_files[1]


contribution_DF <- data.frame()

for (cancer_file in mut_mat_files){
  
  cancer <- gsub("_..*$","",cancer_file)
  
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output/NEW_COHORT")
 
  cancer_mut_mat <- read.table(cancer_file,header=TRUE)
  colnames(cancer_mut_mat) <- gsub("\\.","-",colnames(cancer_mut_mat))
  
  #Normalize each column by its colsums 
  #normalized_mut_mat <- sweep(cancer_mut_mat,2,colSums(cancer_mut_mat),`/`)

  
  #Summarize into one big cohort normalized mutational profile
  #cohort_normalized_mut_mat <- rowSums(normalized_mut_mat)
  
  
  #Alternativley remove high mut samples from samples then rowsums
  
  #test <- colSums(cancer_mut_mat)
  test <- cancer_mut_mat[which(colSums(cancer_mut_mat) >  2000 )]
  #test <- cancer_mut_mat
  test <- rowSums(test)
  test <- as.matrix(test)
  
  fit_to_cosmic <- fit_to_signatures(test,cosmic_signatures)
  
  contribution <- fit_to_cosmic$contribution
  contribution <- contribution[order(contribution,decreasing = TRUE),]
  
  top_x_signatures <- head(contribution,n=10)
  
  append_df <- data.frame("Contribution"=as.vector(top_x_signatures),"Signature"=names(top_x_signatures))
  
  append_df$cancer <- rep(cancer,nrow(append_df))
  
  contribution_DF <- rbind(contribution_DF,append_df)
  
  
  
  

  
}

contribution_DF <- contribution_DF %>% filter(Contribution > 0)

test <- dplyr::count(contribution_DF,Signature)
