#Example pipeline
rm(list=ls())
gc()
#Library
library(nilsPaket)






#Set wd to folder containing all data
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)

#Load files with DNA_repair proteins
DNA_repair <- read.csv("DNA_repair.csv",sep=";")

#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"






#Transform MC3 file to VCF files (needed for mutational_patterns)
#only need to be done once.
if (!dir.exists(file.path(paste(main_wd,"/VCF_files",sep="")))){
  
  dir_name <- paste(main_wd,"/VCF_files",sep="")
  
  dir.create(file.path(dir_name))
  MC3_to_VCF()
  
}  
  


#Read MC3 file
if(!file.exists("MC3.rda")){
  
  MC3 <- read_MC3()
  save(MC3,file="MC3.rda")
  
}
load("MC3.rda")




#Extract highly mutated samples
#Stopped working all of a sudden , cant find Tumor_Sample_Barcode
high_mutation_samples <- get_high_mutations(MC3,"SD",3)


#Subselction of MC3 

#Select genes associated with DNA repair mechanisms
MC3 <- MC3 %>% filter(Hugo_Symbol %in% DNA_repair$hgnc_symbol)

#Select samples with high mutation rates
MC3 <- MC3 %>% filter(Tumor_Sample_Barcode %in% high_mutation_samples)




# Mutational signatures ---------------------------------------------------

vcf_list <- list.files(path=paste(main_wd,"/VCF_files",sep=""),recursive=TRUE,pattern=".vcf")

#select only highly mutated samples
vcf_list_names <- gsub(".vcf","",gsub("[A-Z0-9]*/","",vcf_list))
vcf_list <- vcf_list[vcf_list_names %in% high_mutation_samples]












#Read GISTIC file
file_list <- list.files(pattern="all_data_by_genes.txt",recursive = TRUE)
GISTIC <- read_GISTIC(file_list,DNA_repair$hgnc_symbol)









