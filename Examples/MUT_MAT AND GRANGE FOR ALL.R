#Example pipeline
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


#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
vcf_list_total <- list.files(path=paste(main_wd,"/VCF_files",sep=""),recursive=TRUE,pattern=".vcf")

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output/NEW_cohort")
#Since reading all vcfs at once takes such long time , split it up at cohort level and save to folder


cancer_DF <- data.frame("cancer"=character(0),"number_of_samples"=numeric(0))


#unique(gsub("/..*$","",vcf_list_total))
#"SARC" "SKCM" "STAD" "TGCT" "THCA" "THYM" "UCEC" "UCS"  "UVM"
for (cancer in c("SARC","SKCM" ,"STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS",  "UVM")){
  
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/VCF_files")
  
  print(paste("Starting on ",cancer,sep=""))
  print(Sys.time())
  
  files_to_read <- vcf_list_total[grep(cancer,gsub("/..*$","",vcf_list_total))]
  
  #Remove Cohort VCF-file
  files_to_read <- files_to_read[-grep("/Cohort.vcf",files_to_read)]
  file_names <- gsub(".vcf","",gsub("[A-Z0-9]*/","",files_to_read))
  
  cancer_GRange_VCF <- read_vcfs_as_granges(files_to_read,file_names,ref_genome)
  
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output/NEW_cohort")
  #save(cancer_GRange_VCF,file=paste(cancer,"_GRange.rda",sep=""))
  
  print("Creating mutational matrix")
  
  mutational_matrix_cancer <- mut_matrix(cancer_GRange_VCF,ref_genome)
  #write.table(mutational_matrix_cancer,file=paste(cancer,"_mut_matrix.txt",sep=""))
  
  append_DF <- data.frame("cancer"=cancer,"number_of_samples"=dim(mutational_matrix_cancer)[2])
  cancer_DF <- rbind(cancer_DF,append_DF)
  
}



#Combine all cohort_mut matrices to one big.
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output")
mut_mat_files <- list.files()
mut_mat_files <- mut_mat_files[grep("_mut_matrix.txt",mut_mat_files)]



complete_mut_mat <- do.call(cbind, lapply(mut_mat_files, read.table))
#For some reaseon - turns into . in the TCGA-code
colnames(complete_mut_mat) <- gsub("\\.","-",colnames(complete_mut_mat))


write.table(complete_mut_mat,"COMPLETE_MUT_MAT.txt")


