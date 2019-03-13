#Example pipeline
rm(list=ls())
gc()
#Library
library(nilsPaket)
library(dplyr)

#Set wd to folder containing all data
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)

#Load files with DNA_repair proteins
DNA_repair <- read.csv("DNA_repair.csv",sep=";")



#Read MC3 file
MC3 <- read_MC3()

#Extract highly mutated samples
high_mutation_samples <- get_high_mutations(MC3,type="SD",3)


#Subselction of MC3 

#Select genes associated with DNA repair mechanisms
MC3 <- MC3 %>% filter(Hugo_Symbol %in% DNA_repair$hgnc_symbol)

#Select samples with high mutation rates
MC3 <- MC3 %>% filter(Tumor_Sample_Barcode %in% high_mutation_samples )




#Read GISTIC file
file_list <- list.files(pattern="all_data_by_genes.txt",recursive = TRUE)
GISTIC <- read_GISTIC(file_list,DNA_repair$hgnc_symbol)









