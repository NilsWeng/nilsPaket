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

#Load files with DNA_repair proteins
DNA_repair <- read.csv("DNA_repair.csv",sep=";")

#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"




# MC3 file  ---------------------------------------------------------------------


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
S <- 3
type <- "SD"

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
vcf_list_names <- gsub(".vcf","",gsub("[A-Z0-9]*/","",vcf_list))

#convert VCF files to GRange objects(takes some time)


if (!dir.exists(file.path(paste(main_wd,"/output",sep="")))){
  
  dir_name <- paste(main_wd,"/output",sep="")
  dir.create(file.path(dir_name))
  setwd(paste(main_wd,"/VCF_files",sep=""))
  all_vcfs <- read_vcfs_as_granges(vcf_list,vcf_list_names,ref_genome)
  
  setwd(dir_name)
  save(all_vcfs,file=paste(S,"_all_vcfs.rda",sep=""))
}  

dir_name <- paste(main_wd,"/output",sep="")
load(paste(S,"_all_vcfs.rda",sep=""))



#Create a mutational matrix
setwd(dir_name)
file_name <- paste(type,S,"mut_mat.rda",sep="_")
if(!file.exists(file_name)){
  
  mutational_matrix = mut_matrix(all_vcfs,ref_genome)
  save(mutational_matrix,file=file_name)
} else {
  
  load(file_name)
  
}



#Extract cosine similarity
setwd(main_wd)
#Read cosmic signatures
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))

#Shorten TCGA-barcode
colnames(mutational_matrix) <- gsub("TCGA-","",
                                    gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(mutational_matrix)))

#Create a cosine similarity matrix
cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)
#Cluster cosmic signatures
hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
#store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]
#plot cosine_heatmap
print(plot_cosine_heatmap(cos_sim_samples_cosmic, col_order = cosmic_order, cluster_rows = TRUE))



#Cluster samples based on cosine similarity
N <- 14
sample_cluster <- hclust(dist(cos_sim_samples_cosmic,method="euclidean"),method="complete")
plot_cluster(sample_cluster,N,2.2)





#plot cluster in heatmap
ClusterDF <- plot_cluster_in_cosine(sample_cluster,cos_sim_samples_cosmic,14)

#Manual edit
#Split cluster 3 into 2 
remove <- "13-0889"
new <- c("YC-A89H","DK-A1AC","BT-A2LB")
ClusterDF[ClusterDF$sample %in% new,2] <- "3a"
ClusterDF <- ClusterDF[!ClusterDF$sample %in% remove ,]



#Find signatures in cluster
get_signature(ClusterDF,cos_sim_samples_cosmic,0.6)




#Read GISTIC file
file_list <- list.files(pattern="all_data_by_genes.txt",recursive = TRUE)
GISTIC <- read_GISTIC(file_list,DNA_repair$hgnc_symbol)









