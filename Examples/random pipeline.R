#Example pipeline
rm(list=ls())
gc()
#Library
library(nilsPaket)
library(dplyr)
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(dynamicTreeCut)



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

#Random samples
random_samples <- unique(MC3$Tumor_Sample_Barcode)
random_samples <- sample(random_samples,300)

high_mutation_samples <- random_samples




#Extract highly mutated samples
#Stopped working all of a sudden , cant find Tumor_Sample_Barcode
S <- "random4"
type <- "SD" #SD or absolute

#high_mutation_samples <- get_high_mutations(MC3,type,S)


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




setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Pictures")
pdf(paste(S,".pdf"),height=8.27,width =11.69)





setwd(main_wd)
#Cluster cosmic signatures
hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
#store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]
#plot cosine_heatmap
plot_cosine_heatmap(cos_sim_samples_cosmic, col_order = cosmic_order, cluster_rows = TRUE)



#Cluster samples based on cosine similarity

#N <- 5
sample_cluster <- hclust(dist(cos_sim_samples_cosmic,method="euclidean"),method="complete")

N <- length(unique(cutreeDynamicTree(sample_cluster,maxTreeHeight = 16,deepSplit = TRUE)))

if (N < 5 ){
  
  N <- 10
}


plot_cluster(sample_cluster,N,2.2)



setwd(main_wd)

#plot cluster in heatmap
ClusterDF <- plot_cluster_in_cosine(sample_cluster,cos_sim_samples_cosmic,N)
ClusterDF$sample <- as.character(ClusterDF$sample)
sampleDF$sample <- as.character(sampleDF$sample)
#ClusterDF$sample <- as.numeric(ClusterDF$sample) 
ClusterDF<-left_join(ClusterDF,sampleDF,by="sample")
ClusterDF$TSS.Code <- gsub("TCGA-","",gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",ClusterDF$TCGA_code))
TSS2Study <- read.table("TSS2Studyabb.txt",header=TRUE)
ClusterDF <- left_join(ClusterDF,TSS2Study,by="TSS.Code")
ClusterDF <- ClusterDF %>% select("sample","cluster","TCGA_code","Study.Abbreviation")


load("MC3.rda")
mutDF <- count(MC3,Tumor_Sample_Barcode)
colnames(mutDF) <- c("TCGA_code","mutations")
ClusterDF <- left_join(ClusterDF,mutDF,by="TCGA_code")



rm(TSS2Study,sampleDF,MC3)



#Find signatures in cluster
#Make sure that rownames of cos_sim matches ClusterDF$samples
get_signature(ClusterDF,cos_sim_samples_cosmic,0.6)
library(dplyr)
piechart_cancer_cluster(ClusterDF)




#rm(list=ls()[! ls() %in% c("ClusterDF","DNA_repair","S","type","ref_genome","main_wd")])

# CNV_SNV plot ------------------------------------------------------------


setwd(main_wd)




samples <- ClusterDF %>% select(TCGA_code)
samples <- unique(samples)
samples <- as.character(samples$TCGA_code)
samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)




genes <- read.table("gen_lista.csv",header=TRUE,sep=";")
#genes <- genes %>% filter(System == "BER")
genes <- as.character(genes$Gen)


print(plot_CNV_SNV(samples,genes,S,plot_id = TRUE,cluster_names = TRUE))





####################################################



#Function- Given a set of samples and genes get expression for these


# mRNA Exp ----------------------------------------------------------------

mRNA_DF <- mRNA_exp(genes)
plot_mRNA_SNV(samples,genes,cluster_names=TRUE,mRNA_DF)


setwd(main_wd)
#samples <- ClusterDF %>% filter(cluster %in% c(1,4,8)) %>% select(TCGA_code)

#Select samples



exp <- get_exp_as_matrix(samples,genes,mRNA_DF)
exp <- t(exp)


#For getting cluster lines
test <- ClusterDF 
test$TCGA_code <-  gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",ClusterDF$TCGA_code)
test <- test[which(test$TCGA_code %in% colnames(exp)) ,]
hline_pos <- c()
for (cluster1 in unique(test$cluster)){
  
 
  pos <- tail(which(test$cluster == cluster1),n=1)
  pos <- pos + 0.5
  
  hline_pos <- c(hline_pos,pos)
  
  
}
hline_pos <-head(hline_pos,-1)


#In order to not get grey values
exp[exp >2] <- 2


#Optional: Change TCGA-code to id
test <- test[order(match(test$TCGA_code,colnames(exp))) ,]

colnames(exp) <- paste(rep("_",length(test$sample)),test$sample,sep="")



library(ggplot2)
library(reshape2)

exp_melt <- melt(exp)

p <- ggplot(data = exp_melt, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + labs(x="",y="",title=S)+
  #scale_fill_gradient(low="green", mid="lightblue", high="red",limits
  scale_fill_gradientn(colours=c("red","lightblue","green"), limits=c(0,2))+
  geom_hline(yintercept=hline_pos)+
  
  theme(
    axis.text.x=element_text(colour="black",angle=90, hjust=1,vjust = 0.5,size=6),
    axis.text.y=element_text(vjust = 0.5,size=5)
  )

print(p)

dev.off()
dev.off()


commonly_mut_genes(ClusterDF)






