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
DNA_repair <- read.csv("DNA_repair_extended.csv",sep=";")

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
type <- "SD" #SD or absolute

high_mutation_samples <- get_high_mutations(MC3,type,S)


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

file_name <- paste(S,"_all_vcfs.rda",sep="")
dir_name <- paste(main_wd,"/output",sep="")
setwd(dir_name)

if(!file.exists(file_name)){
  setwd(paste(main_wd,"/VCF_files",sep=""))
  all_vcfs <- read_vcfs_as_granges(vcf_list,vcf_list_names,ref_genome)
  setwd(dir_name)
  save(all_vcfs,file=paste(S,"_all_vcfs.rda",sep=""))
  
} else {
  setwd(dir_name)
  load(file_name)
  
}

#dir_name <- paste(main_wd,"/output",sep="")
#load(paste(S,"_all_vcfs.rda",sep=""))



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


#Just plot number in heatmap.
sampleDF <- data.frame("TCGA_code"=colnames(mutational_matrix))
sampleDF$sample <- c(1:ncol(mutational_matrix))
colnames(mutational_matrix) <- c(1:ncol(mutational_matrix))


#Shorten TCGA-barcode
#colnames(mutational_matrix) <- gsub("TCGA-","",
#                                    gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(mutational_matrix)))

#Create a cosine similarity matrix
cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)

#Cluster cosmic signatures
hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
#store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]
#plot cosine_heatmap
print(plot_cosine_heatmap(cos_sim_samples_cosmic, col_order = cosmic_order, cluster_rows = TRUE))



#Cluster samples based on cosine similarity
N <- 10
sample_cluster <- hclust(dist(cos_sim_samples_cosmic,method="euclidean"),method="complete")
plot_cluster(sample_cluster,N,2.2)





#plot cluster in heatmap
ClusterDF <- plot_cluster_in_cosine(sample_cluster,cos_sim_samples_cosmic,14)
ClusterDF$sample <- as.character(ClusterDF$sample)
sampleDF$sample <- as.character(sampleDF$sample)
#ClusterDF$sample <- as.numeric(ClusterDF$sample) 
ClusterDF<-left_join(ClusterDF,sampleDF,by="sample")
ClusterDF$TSS.Code <- gsub("TCGA-","",gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",ClusterDF$TCGA_code))
TSS2Study <- read.table("TSS2Studyabb.txt",header=TRUE)
ClusterDF <- left_join(ClusterDF,TSS2Study,by="TSS.Code")
ClusterDF <- ClusterDF %>% select("sample","cluster","TCGA_code","Study.Abbreviation")


load("MC3.rda")
mutDF <- count(MC3,"Tumor_Sample_Barcode")
colnames(mutDF) <- c("TCGA_code","mutations")
ClusterDF <- left_join(ClusterDF,mutDF,by="TCGA_code")



rm(TSS2Study,sampleDF,MC3)


#Manual edit
#Split cluster 3 into 2 
remove <- "TCGA-13-0889-01A-01W-0420-08"
new <- c("TCGA-YC-A89H-01A-11D-A364-08","TCGA-DK-A1AC-01A-11D-A13W-08","TCGA-BT-A2LB-01A-11D-A18F-08")
ClusterDF[ClusterDF$sample %in% new,2] <- "3a"
ClusterDF <- ClusterDF[!ClusterDF$TCGA_code %in% remove ,]


#Find signatures in cluster
#Make sure that rownames of cos_sim matches ClusterDF$samples
get_signature(ClusterDF,cos_sim_samples_cosmic,0.6)

rm(list=ls()[! ls() %in% c("ClusterDF","DNA_repair","S","type","ref_genome","main_wd")])

# CNV_SNV plot ------------------------------------------------------------


setwd(main_wd)

samples <- ClusterDF %>% filter(cluster == 3) %>% select(TCGA_code)
samples <- as.character(samples$TCGA_code)
samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)


genes <- as.character(DNA_repair$hgnc_symbol)

#################################################################################
#Test MMR cluster 6,7,9
samples <- ClusterDF %>% filter(cluster %in% c(2)) %>% select(TCGA_code)
samples <- as.character(samples$TCGA_code)
samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)



genes <- read.table("gen_lista.csv",header=TRUE,sep=";")
genes <- genes %>% filter(System == "NER")
genes <- as.character(genes$Gen)

plot_CNV_SNV(samples,genes,"TEST")



#Just take random genes 
samples <- as.character(sample(MC3$Tumor_Sample_Barcode,40))
samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)
genes <- as.character(unique(MC3$Hugo_Symbol[1:40]))
plot_CNV_SNV(samples,genes,"Random")

####################################################


#Plot and save all cluster
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Pictures/CNV_SNV")
pdf(file="one_cluster_a_time.pdf",width=16, height=10)

for (Cluster in unique(ClusterDF$cluster)){
  
  setwd(main_wd)
  
  samples <- ClusterDF %>% filter(cluster == Cluster) %>% select(TCGA_code)
  samples <- as.character(samples$TCGA_code)
  samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)
  
  
  genes <- as.character(DNA_repair$hgnc_symbol)
  print(plot_CNV_SNV(samples,genes,paste("Cluster",Cluster,sep=" ")))
  
  
}

dev.off()



#Function- Given a set of samples and genes get expression for these


mRNA_exp <- function(samples,genes){
  
  setwd(main_wd)
  file_list <- list.files(pattern="_normalized__data.data.txt",recursive = TRUE)
  
  read_mRNASeq <- function(files){
    
    
    
    Table <- fread(files,header=TRUE,stringsAsFactors = FALSE)
    #Filter only for genes involved in DNA-repair
    Table$`Hybridization REF` <- gsub("\\|.*$","",Table$`Hybridization REF`)
    #Table <- Table %>% filter(Table$`Hybridization REF` %in% genes)
    
    print("done with")
    return(Table)
    
    
  }
  
  #file_list <- file_list[1:3]
  
  dataset <- do.call("cbind",lapply(file_list ,read_mRNASeq))
  dataset <- dataset[, !duplicated(colnames(dataset))]
  
  setwd(main_wd)
  write.table(dataset,"PANCAN_mRNA_expression.txt")
  save(dataset,file="PANCAN_mRNA_expression.rda")
  
  
}



get_cancer_abb <- function(TCGA_code){
  
  setwd(main_wd)
  TSS2Study <- read.table("TSS2Studyabb.txt",header=TRUE)
  
  TSS_code <- gsub("TCGA-","",gsub(""))
  
  
  
  
}



