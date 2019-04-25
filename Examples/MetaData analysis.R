#Script for getting Metadata about the different TCGA Studies. I.e how many samples an so on
rm(list=ls())
gc()
#Library
library(nilsPaket)
library(data.table)

main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)




#MC3-file
if(!file.exists("MC3.rda")){
  
  MC3 <- read_MC3()
  save(MC3,file="MC3.rda")
  
}
load("MC3.rda")
MC3_samples <- unique(MC3$Tumor_Sample_Barcode)
rm(MC3)




#GISTIC files
gistic_files <- list.files(pattern="all_thresholded.by_genes",recursive = TRUE)
gistic_samples <- c()
for (file_name in gistic_files){
  
  gistic_file <- fread(file_name,header=TRUE,stringsAsFactors = FALSE)
  gistic_samples_cancer <- colnames(gistic_file)
  gistic_samples_cancer <- gistic_samples_cancer[4:length(gistic_samples_cancer)]

  gistic_samples <- c(gistic_samples,gistic_samples_cancer)
  
}
rm(gistic_file,file_name,gistic_samples_cancer,gistic_files)



#mRNA-expression files

mRNA_files <- list.files(pattern="RSEM_genes_normalized__data.data",recursive = TRUE)

mRNA_samples <- c()
for (file_name in mRNA_files){
  
  mRNA_file <- fread(file_name,header=TRUE,stringsAsFactors = FALSE)
  
  mRNA_samples_cancer <- colnames(mRNA_file)
  mRNA_samples_cancer <- mRNA_samples_cancer[2:length(mRNA_samples_cancer)]
  
  mRNA_samples <- c(mRNA_samples,mRNA_samples_cancer)
  
}

rm(mRNA_file,file_name,mRNA_samples_cancer,mRNA_files)



#Compare studies
#Normalize TCGA-barcode
gistic_samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",gistic_samples)
MC3_samples    <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",MC3_samples)
mRNA_samples   <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",mRNA_samples)


#Some samples in mRNA are normal samples , remove these
sample_type <- gsub("[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-","",mRNA_samples)
sample_type <- as.numeric(gsub("[A-Z]$","",sample_type))
delete_sampels <- which(sample_type > 9)
mRNA_samples <- mRNA_samples[-delete_sampels]





shared <- Reduce(intersect, list(MC3_samples,mRNA_samples,gistic_samples))
shared <- get_cancer_type(shared)
shared <- as.data.frame(table(shared))




gistic_cancer <- get_cancer_type(gistic_samples)
gistic_cancer <- as.data.frame(table(gistic_cancer))


MC3_cancer <- get_cancer_type(MC3_samples)
MC3_cancer <- as.data.frame(table(MC3_cancer))


mRNA_cancer <- get_cancer_type(mRNA_samples)
mRNA_cancer <- as.data.frame(table(mRNA_cancer))



 




metadata_DF <- data.frame("Study_Abbreviation" = MC3_cancer$MC3_cancer,"MC3" = MC3_cancer$Freq,"GISTIC" = gistic_cancer$Freq,"mRNA-expression" = mRNA_cancer$Freq,"Shared" = shared$Freq)


StudyName <- read.table("StudyName.csv",header=TRUE,sep=";")
colnames(StudyName) <- c("Study_Abbreviation","Study Name")

library(dplyr)
metadata_DF <- left_join(StudyName,metadata_DF)


tot_vect <- c("Total","-",sum(metadata_DF$MC3),sum(metadata_DF$GISTIC),
              sum(metadata_DF$`mRNA expression`),sum(metadata_DF$`Shared between data types`))

tot_vect <- as.data.frame(tot_vect)
tot_vect <- as.data.frame(t(tot_vect))
colnames(tot_vect) <- colnames(metadata_DF)

metadata_DF <- rbind(metadata_DF,tot_vect)




colnames(metadata_DF) <- c("Study Abrbeviation","Study Name","MC3","GISTIC","mRNA expression","Shared between all data types")
library(formattable)
print(formattable(metadata_DF,align="l"))


library(gridExtra)
library(grid)


pdf(file = "myfile.pdf", height = 12, width = 26)
grid.newpage()
grid.table(metadata_DF[1:34, ],rows = NULL)
dev.off()

